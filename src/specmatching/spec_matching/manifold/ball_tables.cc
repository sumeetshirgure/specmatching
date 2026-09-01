// Copyright 2026 SpecMatching contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "specmatching/spec_matching/manifold/ball_tables.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstring>
#include <limits>
#include <queue>
#include <stdexcept>
#include <thread>

namespace pm {
namespace spec_matching {

namespace {

constexpr uint32_t NO_PRED = std::numeric_limits<uint32_t>::max();
constexpr pm::cumulative_time_int INF_DIST = std::numeric_limits<pm::cumulative_time_int>::max();

/// One source node's slice of the tables, in local (per-node) CSR form. Kept flat rather than as a
/// vector of vectors: at `d = 21` a ball holds thousands of entries and a heap allocation per entry
/// would dominate the compile.
struct NodeBall {
    std::vector<uint32_t> target;
    std::vector<pm::weight_int> w_int;
    std::vector<uint64_t> mask_offsets{0};
    std::vector<uint32_t> mask_ids;
    std::vector<uint64_t> path_offsets{0};
    std::vector<uint32_t> path_nodes;

    bool has_bcost{false};
    pm::weight_int bcost{0};
    std::vector<uint32_t> bcost_mask;
    std::vector<uint32_t> bcost_path;

    /// Mirrors `BallParams::store_paths`; kept per-ball so the worker does not have to reach back
    /// into the params on the inner loop.
    bool store_paths{false};
    uint64_t ambiguous{0};
};

/// Per-thread Dijkstra state, reset through a touched list so that a sweep over all sources costs
/// nothing proportional to the graph per source.
struct DijkstraScratch {
    std::vector<pm::cumulative_time_int> dist;
    std::vector<uint32_t> pred;
    std::vector<uint8_t> settled;
    std::vector<uint32_t> touched;
    std::vector<uint32_t> reached;
    std::priority_queue<
        std::pair<pm::cumulative_time_int, uint32_t>,
        std::vector<std::pair<pm::cumulative_time_int, uint32_t>>,
        std::greater<std::pair<pm::cumulative_time_int, uint32_t>>>
        frontier;

    explicit DijkstraScratch(size_t num_nodes)
        : dist(num_nodes, INF_DIST), pred(num_nodes, NO_PRED), settled(num_nodes, 0) {
    }

    void reset() {
        for (uint32_t v : touched) {
            dist[v] = INF_DIST;
            pred[v] = NO_PRED;
            settled[v] = 0;
        }
        touched.clear();
        reached.clear();
        while (!frontier.empty())
            frontier.pop();
    }
};

/// XORs a sorted observable id list into a sorted accumulator.
void xor_merge_ids(std::vector<uint32_t>& accumulator, const std::vector<size_t>& ids, std::vector<uint32_t>& scratch) {
    scratch.clear();
    size_t i = 0;
    size_t j = 0;
    while (i < accumulator.size() && j < ids.size()) {
        uint32_t a = accumulator[i];
        uint32_t b = (uint32_t)ids[j];
        if (a < b) {
            scratch.push_back(a);
            i++;
        } else if (b < a) {
            scratch.push_back(b);
            j++;
        } else {
            // Present in both: crossed twice, so it cancels.
            i++;
            j++;
        }
    }
    while (i < accumulator.size())
        scratch.push_back(accumulator[i++]);
    while (j < ids.size())
        scratch.push_back((uint32_t)ids[j++]);
    accumulator.swap(scratch);
}

/// Everything the per-source sweep needs, shared read-only across threads.
struct BallCompiler {
    const pm::MatchingGraph& graph;
    const pm::SearchGraph* search_graph;
    /// Per node, the indices into `neighbors` of its non-boundary neighbours, sorted by target node
    /// id. Relaxing in this order is half of the canonical-path rule; the other half is the
    /// smaller-predecessor tie-break below.
    std::vector<std::vector<uint32_t>> sorted_neighbors;
    horizon_int r_int;
    bool use_obs_masks;

    BallCompiler(const pm::MatchingGraph& graph, const pm::SearchGraph* search_graph, horizon_int r_int)
        : graph(graph),
          search_graph(search_graph),
          r_int(r_int),
          use_obs_masks(graph.num_observables <= sizeof(pm::obs_int) * 8) {
        const pm::DetectorNode* base = graph.nodes.data();
        sorted_neighbors.resize(graph.nodes.size());
        for (size_t v = 0; v < graph.nodes.size(); v++) {
            const auto& node = graph.nodes[v];
            auto& order = sorted_neighbors[v];
            for (size_t i = 0; i < node.neighbors.size(); i++) {
                if (node.neighbors[i] == nullptr)
                    continue;  // the boundary half-edge, handled separately
                order.push_back((uint32_t)i);
            }
            std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
                return node.neighbors[a] - base < node.neighbors[b] - base;
            });
        }
    }

    inline uint32_t neighbor_id(size_t v, uint32_t index) const {
        return (uint32_t)(graph.nodes[v].neighbors[index] - graph.nodes.data());
    }

    /// Bounded Dijkstra from `src` on the **integer** weights. Never on a float distance: sparse
    /// blossom sums rounded edge weights, so `round(sum of floats)` would be a different metric
    /// (§M2.4), and the identity test would fail on a minority of shots with no pattern.
    ///
    /// `prefer_larger_pred` inverts the canonical tie-break; that run is the mask certification,
    /// never the stored table.
    void run_dijkstra(DijkstraScratch& s, uint32_t src, bool prefer_larger_pred) const {
        s.reset();
        s.dist[src] = 0;
        s.touched.push_back(src);
        s.frontier.push({0, src});
        while (!s.frontier.empty()) {
            auto [d, u] = s.frontier.top();
            s.frontier.pop();
            if (s.settled[u] || d != s.dist[u])
                continue;
            s.settled[u] = 1;
            if (u != src)
                s.reached.push_back(u);

            const auto& order = sorted_neighbors[u];
            size_t degree = order.size();
            for (size_t k = 0; k < degree; k++) {
                uint32_t index = order[prefer_larger_pred ? degree - 1 - k : k];
                uint32_t v = neighbor_id(u, index);
                pm::cumulative_time_int candidate = d + (pm::cumulative_time_int)graph.nodes[u].neighbor_weights[index];
                if (candidate > r_int)
                    continue;
                if (candidate < s.dist[v]) {
                    if (s.dist[v] == INF_DIST)
                        s.touched.push_back(v);
                    s.dist[v] = candidate;
                    s.pred[v] = u;
                    s.frontier.push({candidate, v});
                } else if (candidate == s.dist[v] && v != src) {
                    // Every valid predecessor of `v` is settled no later than `v` itself, so by the
                    // end of the sweep `pred[v]` is the extremal one over *all* of them. That is
                    // what makes the path canonical rather than dependent on the queue order.
                    bool better = prefer_larger_pred ? (s.pred[v] != NO_PRED && u > s.pred[v]) : (u < s.pred[v]);
                    if (better)
                        s.pred[v] = u;
                }
            }
        }
    }

    /// The observable ids crossed by the canonical path from `src` to `u`, and (optionally) the
    /// node sequence itself. `path_out` is written from `src` to `u` inclusive.
    void reconstruct(
        const DijkstraScratch& s,
        uint32_t src,
        uint32_t u,
        std::vector<uint32_t>& mask_out,
        std::vector<uint32_t>* path_out,
        std::vector<uint32_t>& scratch) const {
        mask_out.clear();
        if (path_out != nullptr)
            path_out->clear();

        pm::obs_int mask = 0;
        uint32_t node = u;
        while (node != src) {
            uint32_t p = s.pred[node];
            assert(p != NO_PRED && "canonical path ran off the predecessor tree");
            if (path_out != nullptr)
                path_out->push_back(node);
            size_t index = graph.nodes[p].index_of_neighbor(const_cast<pm::DetectorNode*>(&graph.nodes[node]));
            if (use_obs_masks) {
                mask ^= graph.nodes[p].neighbor_observables[index];
            } else {
                const auto& search_node = search_graph->nodes[p];
                size_t search_index =
                    search_node.index_of_neighbor(const_cast<pm::SearchDetectorNode*>(&search_graph->nodes[node]));
                xor_merge_ids(mask_out, search_node.neighbor_observable_indices[search_index], scratch);
            }
            node = p;
        }
        if (path_out != nullptr) {
            path_out->push_back(src);
            std::reverse(path_out->begin(), path_out->end());
        }
        if (use_obs_masks) {
            for (uint32_t bit = 0; mask != 0; bit++, mask >>= 1) {
                if (mask & 1)
                    mask_out.push_back(bit);
            }
        }
    }

    /// The cheapest boundary path from `src`, over the nodes the bounded sweep settled. Exact
    /// whenever the answer is `<= R`: any boundary path of cost `<= R` has every prefix `<= R`, so
    /// its last node was settled.
    bool find_boundary(
        const DijkstraScratch& s, uint32_t src, uint32_t& node_out, pm::cumulative_time_int& cost_out) const {
        pm::cumulative_time_int best = INF_DIST;
        uint32_t best_node = NO_PRED;
        auto consider = [&](uint32_t v) {
            const auto& node = graph.nodes[v];
            if (node.neighbors.empty() || node.neighbors[0] != nullptr)
                return;
            pm::cumulative_time_int candidate = s.dist[v] + (pm::cumulative_time_int)node.neighbor_weights[0];
            if (candidate > r_int)
                return;
            if (candidate < best || (candidate == best && v < best_node)) {
                best = candidate;
                best_node = v;
            }
        };
        consider(src);
        for (uint32_t v : s.reached)
            consider(v);
        if (best_node == NO_PRED)
            return false;
        node_out = best_node;
        cost_out = best;
        return true;
    }

    void compile_node(DijkstraScratch& s, uint32_t src, NodeBall& out, bool certify) const {
        run_dijkstra(s, src, false);

        std::vector<uint32_t> entries = s.reached;
        std::sort(entries.begin(), entries.end(), [&](uint32_t a, uint32_t b) {
            if (s.dist[a] != s.dist[b])
                return s.dist[a] < s.dist[b];
            return a < b;
        });

        std::vector<uint32_t> mask;
        std::vector<uint32_t> path;
        std::vector<uint32_t> scratch;

        out.target.reserve(entries.size());
        out.w_int.reserve(entries.size());
        for (uint32_t u : entries) {
            out.target.push_back(u);
            out.w_int.push_back((pm::weight_int)s.dist[u]);
            reconstruct(s, src, u, mask, out.store_paths ? &path : nullptr, scratch);
            out.mask_ids.insert(out.mask_ids.end(), mask.begin(), mask.end());
            out.mask_offsets.push_back(out.mask_ids.size());
            if (out.store_paths) {
                out.path_nodes.insert(out.path_nodes.end(), path.begin(), path.end());
                out.path_offsets.push_back(out.path_nodes.size());
            }
        }

        uint32_t boundary_node = 0;
        pm::cumulative_time_int boundary_cost = 0;
        if (find_boundary(s, src, boundary_node, boundary_cost)) {
            out.has_bcost = true;
            out.bcost = (pm::weight_int)boundary_cost;
            reconstruct(s, src, boundary_node, mask, out.store_paths ? &path : nullptr, scratch);
            // The half-edge into the boundary itself is part of the path's homology.
            pm::obs_int half_edge_mask = graph.nodes[boundary_node].neighbor_observables[0];
            if (use_obs_masks) {
                std::vector<uint32_t> half_edge_ids;
                for (uint32_t bit = 0; half_edge_mask != 0; bit++, half_edge_mask >>= 1) {
                    if (half_edge_mask & 1)
                        half_edge_ids.push_back(bit);
                }
                std::vector<size_t> as_size_t(half_edge_ids.begin(), half_edge_ids.end());
                xor_merge_ids(mask, as_size_t, scratch);
            } else {
                const auto& search_node = search_graph->nodes[boundary_node];
                size_t half_edge_index = search_node.index_of_neighbor(nullptr);
                xor_merge_ids(mask, search_node.neighbor_observable_indices[half_edge_index], scratch);
            }
            out.bcost_mask = mask;
            if (out.store_paths)
                out.bcost_path = path;
        }

        if (!certify)
            return;

        // Certification: the same tables under a perturbed tie-break. A differing observable set
        // means the pair has two shortest paths of differing homology — recorded, never blocking.
        std::vector<uint32_t> alt_mask;
        run_dijkstra(s, src, true);
        for (size_t e = 0; e < out.target.size(); e++) {
            uint32_t u = out.target[e];
            assert(s.dist[u] == (pm::cumulative_time_int)out.w_int[e] && "perturbed tie-break changed a distance");
            reconstruct(s, src, u, alt_mask, nullptr, scratch);
            size_t begin = out.mask_offsets[e];
            size_t end = out.mask_offsets[e + 1];
            bool same = (end - begin) == alt_mask.size() &&
                        std::equal(out.mask_ids.begin() + begin, out.mask_ids.begin() + end, alt_mask.begin());
            if (!same)
                out.ambiguous++;
        }
    }
};

}  // namespace

uint64_t hash_matching_graph(const pm::MatchingGraph& graph) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](uint64_t value) {
        for (int byte = 0; byte < 8; byte++) {
            h ^= (value >> (8 * byte)) & 0xff;
            h *= 1099511628211ull;
        }
    };
    mix(graph.nodes.size());
    mix(graph.num_observables);
    uint64_t normalising_bits;
    static_assert(sizeof(normalising_bits) == sizeof(graph.normalising_constant), "double is not 8 bytes");
    std::memcpy(&normalising_bits, &graph.normalising_constant, sizeof(normalising_bits));
    mix(normalising_bits);
    const pm::DetectorNode* base = graph.nodes.data();
    for (const auto& node : graph.nodes) {
        mix(node.neighbors.size());
        for (size_t i = 0; i < node.neighbors.size(); i++) {
            mix(node.neighbors[i] == nullptr ? UINT64_MAX : (uint64_t)(node.neighbors[i] - base));
            mix(node.neighbor_weights[i]);
            mix(node.neighbor_observables[i]);
        }
    }
    for (bool is_boundary : graph.is_user_graph_boundary_node)
        mix(is_boundary ? 1 : 0);
    return h;
}

BallTables compile_ball_tables(const pm::Mwpm& mwpm, const BallParams& params, size_t num_threads) {
    params.validate();
    auto wall_start = std::chrono::high_resolution_clock::now();

    const pm::MatchingGraph& graph = mwpm.flooder.graph;
    const pm::SearchGraph* search_graph = nullptr;
    if (graph.num_observables > sizeof(pm::obs_int) * 8) {
        if (mwpm.search_flooder.graph.nodes.size() != graph.nodes.size())
            throw std::invalid_argument(
                "compile_ball_tables needs a search graph when the DEM has more than 64 observables: the obs_int "
                "masks the matching graph stores are only filled below that threshold.");
        search_graph = &mwpm.search_flooder.graph;
    }

    BallTables tables;
    tables.params = params;
    tables.normalising_constant = graph.normalising_constant;
    tables.t_max_int = to_time_units(params.T_max, graph.normalising_constant);
    // Invariant 6 is a statement in weight units; the time-unit form is derived. Deriving it by two
    // *independent* roundings does not preserve it — `round(2x)` can be `2 * round(x) - 1` — so at
    // the exact-equality setting `R = 2 * T_max`, which is the natural thing for a caller to write
    // and what the python helper defaults to, the derived form can fail by one time unit.
    //
    // Widening the radius is the safe direction: a ball that is one unit too big costs memory,
    // whereas one that is one unit too small silently omits a reachable pair and changes the answer
    // (§M2.0). `BallParams::validate` has already rejected `R < 2 * T_max` in weight units, which is
    // the check that actually protects the theorem; this only stops a rounding artifact from
    // presenting as a configuration error.
    tables.r_int = std::max(to_time_units(params.R, graph.normalising_constant), 2 * tables.t_max_int);
    tables.num_nodes = graph.nodes.size();
    tables.num_observables = graph.num_observables;
    tables.graph_hash = hash_matching_graph(graph);

    tables.shell_width_int = params.shell_width > 0 ? to_time_units(params.shell_width, graph.normalising_constant) : 0;
    if (params.shell_width > 0 && tables.shell_width_int <= 0)
        throw std::invalid_argument("BallParams::shell_width is below one time unit of the discretised metric.");
    tables.num_shells = tables.shell_width_int > 0 ? (size_t)(tables.r_int / tables.shell_width_int) + 1 : 1;

    size_t num_nodes = graph.nodes.size();
    BallCompiler compiler(graph, search_graph, tables.r_int);

    std::vector<NodeBall> per_node(num_nodes);
    for (auto& ball : per_node)
        ball.store_paths = params.store_paths;

    if (num_threads == 0)
        num_threads = std::max<size_t>(1, std::thread::hardware_concurrency());
    num_threads = std::min(num_threads, std::max<size_t>(1, num_nodes));

    // Embarrassingly parallel over source nodes. The output is independent of the thread count:
    // each source's slice is written to its own `NodeBall` and the pools are concatenated in
    // ascending node id below, so compilation stays bit-reproducible (invariant: determinism, §0).
    auto worker = [&](size_t thread_index) {
        DijkstraScratch scratch(num_nodes);
        for (size_t v = thread_index; v < num_nodes; v += num_threads)
            compiler.compile_node(scratch, (uint32_t)v, per_node[v], params.certify_masks);
    };
    if (num_threads == 1) {
        worker(0);
    } else {
        std::vector<std::thread> threads;
        threads.reserve(num_threads);
        for (size_t t = 0; t < num_threads; t++)
            threads.emplace_back(worker, t);
        for (auto& thread : threads)
            thread.join();
    }

    // ---- Concatenate into the CSR pools, in ascending node id.
    tables.ball_offsets.reserve(num_nodes + 1);
    tables.ball_offsets.push_back(0);
    tables.ball_mask_offsets.push_back(0);
    if (params.store_paths)
        tables.ball_path_offsets.push_back(0);
    tables.has_bcost.resize(num_nodes, 0);
    tables.bcost_w_int.resize(num_nodes, 0);
    tables.bcost_mask_offsets.reserve(num_nodes + 1);
    tables.bcost_mask_offsets.push_back(0);
    if (params.store_paths) {
        tables.bcost_path_offsets.reserve(num_nodes + 1);
        tables.bcost_path_offsets.push_back(0);
    }
    tables.ball_shell_offsets.resize(num_nodes * (tables.num_shells + 1));
    tables.stats.shell_entry_counts.assign(tables.num_shells, 0);

    for (size_t v = 0; v < num_nodes; v++) {
        NodeBall& ball = per_node[v];
        uint64_t base = tables.ball_target.size();
        tables.ball_target.insert(tables.ball_target.end(), ball.target.begin(), ball.target.end());
        tables.ball_w_int.insert(tables.ball_w_int.end(), ball.w_int.begin(), ball.w_int.end());
        uint64_t mask_base = tables.ball_mask_ids.size();
        tables.ball_mask_ids.insert(tables.ball_mask_ids.end(), ball.mask_ids.begin(), ball.mask_ids.end());
        for (size_t e = 1; e < ball.mask_offsets.size(); e++)
            tables.ball_mask_offsets.push_back(mask_base + ball.mask_offsets[e]);
        if (params.store_paths) {
            uint64_t path_base = tables.ball_path_nodes.size();
            tables.ball_path_nodes.insert(tables.ball_path_nodes.end(), ball.path_nodes.begin(), ball.path_nodes.end());
            for (size_t e = 1; e < ball.path_offsets.size(); e++)
                tables.ball_path_offsets.push_back(path_base + ball.path_offsets[e]);
        }
        tables.ball_offsets.push_back(tables.ball_target.size());

        // Shell boundaries. Entries are sorted by `(w_int, target)`, so each shell is a contiguous
        // run and the boundaries are found by a single walk.
        size_t stride = tables.num_shells + 1;
        size_t cursor = 0;
        for (size_t s = 0; s <= tables.num_shells; s++) {
            if (s == tables.num_shells) {
                tables.ball_shell_offsets[v * stride + s] = base + ball.target.size();
                break;
            }
            tables.ball_shell_offsets[v * stride + s] = base + cursor;
            pm::cumulative_time_int shell_end = tables.shell_width_int > 0
                                                    ? (pm::cumulative_time_int)(s + 1) * tables.shell_width_int
                                                    : tables.r_int + 1;
            size_t run_start = cursor;
            while (cursor < ball.w_int.size() && (pm::cumulative_time_int)ball.w_int[cursor] < shell_end)
                cursor++;
            tables.stats.shell_entry_counts[s] += cursor - run_start;
        }

        if (ball.has_bcost) {
            tables.has_bcost[v] = 1;
            tables.bcost_w_int[v] = ball.bcost;
            tables.stats.nodes_with_boundary++;
        }
        tables.bcost_mask_ids.insert(tables.bcost_mask_ids.end(), ball.bcost_mask.begin(), ball.bcost_mask.end());
        tables.bcost_mask_offsets.push_back(tables.bcost_mask_ids.size());
        if (params.store_paths) {
            tables.bcost_path_nodes.insert(
                tables.bcost_path_nodes.end(), ball.bcost_path.begin(), ball.bcost_path.end());
            tables.bcost_path_offsets.push_back(tables.bcost_path_nodes.size());
        }

        tables.stats.ambiguous_mask_pairs += ball.ambiguous;
        tables.stats.max_ball_size = std::max(tables.stats.max_ball_size, (uint64_t)ball.target.size());

        // Free the per-node slice as we go: the concatenated pools are already the peak, and
        // holding both costs twice the memory at `d = 21`.
        ball = NodeBall();
    }

    // ---- Bitset view.
    tables.ball_word_offsets.reserve(num_nodes + 1);
    tables.ball_word_offsets.push_back(0);
    tables.ball_word_base.resize(num_nodes, 0);
    tables.ball_entry_by_rank.resize(tables.ball_target.size(), 0);
    std::vector<std::pair<uint32_t, uint64_t>> by_target;
    for (size_t v = 0; v < num_nodes; v++) {
        uint64_t begin = tables.ball_offsets[v];
        uint64_t end = tables.ball_offsets[v + 1];
        if (begin == end) {
            tables.ball_word_offsets.push_back(tables.ball_words.size());
            continue;
        }
        by_target.clear();
        by_target.reserve(end - begin);
        for (uint64_t e = begin; e < end; e++)
            by_target.emplace_back(tables.ball_target[e], e);
        std::sort(by_target.begin(), by_target.end());
        for (size_t k = 0; k < by_target.size(); k++)
            tables.ball_entry_by_rank[begin + k] = by_target[k].second;

        uint32_t word_base = by_target.front().first / 64;
        uint32_t word_last = by_target.back().first / 64;
        tables.ball_word_base[v] = word_base;
        size_t words_start = tables.ball_words.size();
        tables.ball_words.resize(words_start + (word_last - word_base + 1), 0);
        for (const auto& [target, weight] : by_target)
            tables.ball_words[words_start + target / 64 - word_base] |= (uint64_t)1 << (target % 64);
        tables.ball_word_rank.resize(tables.ball_words.size(), 0);
        uint32_t running = 0;
        for (size_t k = words_start; k < tables.ball_words.size(); k++) {
            tables.ball_word_rank[k] = running;
            running += (uint32_t)__builtin_popcountll(tables.ball_words[k]);
        }
        assert(running == by_target.size() && "bitset view lost an entry");
        tables.ball_word_offsets.push_back(tables.ball_words.size());
    }

    // ---- Stats.
    auto& stats = tables.stats;
    stats.num_nodes = num_nodes;
    stats.total_entries = tables.ball_target.size();
    stats.mean_ball_size = num_nodes ? (double)stats.total_entries / (double)num_nodes : 0;
    uint64_t total_words = 0;
    for (size_t v = 0; v < num_nodes; v++) {
        uint64_t len = tables.word_len(v);
        total_words += len;
        stats.max_ball_word_len = std::max(stats.max_ball_word_len, len);
    }
    stats.mean_ball_word_len = num_nodes ? (double)total_words / (double)num_nodes : 0;

    auto bytes_of = [](const auto& vec) {
        return (uint64_t)(vec.size() * sizeof(typename std::decay_t<decltype(vec)>::value_type));
    };
    stats.bytes_targets = bytes_of(tables.ball_target);
    stats.bytes_weights = bytes_of(tables.ball_w_int);
    stats.bytes_masks = bytes_of(tables.ball_mask_ids) + bytes_of(tables.ball_mask_offsets) +
                        bytes_of(tables.bcost_mask_ids) + bytes_of(tables.bcost_mask_offsets);
    stats.bytes_paths = bytes_of(tables.ball_path_nodes) + bytes_of(tables.ball_path_offsets) +
                        bytes_of(tables.bcost_path_nodes) + bytes_of(tables.bcost_path_offsets);
    stats.bytes_words = bytes_of(tables.ball_words) + bytes_of(tables.ball_word_offsets) +
                        bytes_of(tables.ball_word_base) + bytes_of(tables.ball_word_rank) +
                        bytes_of(tables.ball_entry_by_rank);
    stats.bytes_boundary = bytes_of(tables.has_bcost) + bytes_of(tables.bcost_w_int);
    stats.bytes_offsets = bytes_of(tables.ball_offsets) + bytes_of(tables.ball_shell_offsets);
    stats.bytes_total = stats.bytes_targets + stats.bytes_weights + stats.bytes_masks + stats.bytes_paths +
                        stats.bytes_words + stats.bytes_boundary + stats.bytes_offsets;
    stats.compile_wall_seconds =
        std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - wall_start).count();

#ifndef NDEBUG
    // Invariant 7: symmetry. `u in B(v) <=> v in B(u)`, with equal `w_int`. A violation means the
    // Dijkstra bound or the tie-break is asymmetric, which would make `H` depend on which endpoint
    // of a pair the shot-time intersection happened to walk from.
    for (size_t v = 0; v < num_nodes; v++) {
        for (uint64_t e = tables.ball_offsets[v]; e < tables.ball_offsets[v + 1]; e++) {
            uint32_t u = tables.ball_target[e];
            uint64_t begin = tables.ball_offsets[u];
            uint64_t end = tables.ball_offsets[u + 1];
            bool found = false;
            for (uint64_t f = begin; f < end; f++) {
                if (tables.ball_target[f] == v) {
                    assert(tables.ball_w_int[f] == tables.ball_w_int[e] && "ball weights are not symmetric");
                    found = true;
                    break;
                }
            }
            assert(found && "ball membership is not symmetric");
        }
    }
#endif

    return tables;
}

}  // namespace spec_matching
}  // namespace pm
