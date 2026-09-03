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

#include "specmatching/spec_matching/manifold/ball_graph.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>

#include "specmatching/spec_matching/perf/spec_matching_profile.h"

namespace pm {
namespace spec_matching {

namespace {

constexpr uint32_t NOT_A_DEFECT = std::numeric_limits<uint32_t>::max();

/// The sum of every buffer's capacity. Comparing it before and after a build is a cheap, allocator
/// -independent way to notice that the arena grew.
uint64_t total_capacity(const BallGraphArena& arena) {
    return arena.graph.h_to_det.capacity() + arena.graph.edges.capacity() + arena.graph.boundary_edges.capacity() +
           arena.det_to_h.capacity() + arena.syndrome_words.capacity() + arena.touched_words.capacity();
}

/// Bytes of the observable id list behind one ball entry, plus the two CSR offsets that bound it.
/// Physically these are read by `BallMwpm::rebuild`, not by the loop below; they are charged to the
/// intersection because the pair is what selects them, and because the hardware stage the counter
/// is sizing does the selection and the fetch together.
uint64_t mask_bytes(const BallTables& tables, uint64_t entry) {
    using namespace ball_element_bytes;
    return 2 * MASK_OFFSET + (tables.ball_mask_offsets[entry + 1] - tables.ball_mask_offsets[entry]) * MASK_ID;
}

uint64_t boundary_mask_bytes(const BallTables& tables, uint64_t det) {
    using namespace ball_element_bytes;
    return 2 * BCOST_MASK_OFFSET +
           (tables.bcost_mask_offsets[det + 1] - tables.bcost_mask_offsets[det]) * BCOST_MASK_ID;
}

/// What `BITSET` would have read for `det`, without running it: the node's bitset window, clamped
/// exactly where the loop clamps it — a window whose tail runs past the syndrome bitset is cut
/// short there.
uint64_t bitset_words_counterfactual(const BallTables& tables, uint64_t det, size_t syndrome_words) {
    uint64_t base = tables.ball_word_base[det];
    if (base >= syndrome_words)
        return 0;
    return std::min(tables.word_len(det), (uint64_t)syndrome_words - base);
}

/// What `SCAN` would have walked for `det` at this horizon, without running it. The entries are
/// sorted by `(w_int, target)`, so the walk stops at the first one past `2T` — and that entry's
/// weight is still read, which is the `+ 1`.
uint64_t scan_entries_counterfactual(const BallTables& tables, uint64_t det, pm::cumulative_time_int two_t) {
    auto first = tables.ball_w_int.begin() + (ptrdiff_t)tables.ball_begin(det);
    auto last = tables.ball_w_int.begin() + (ptrdiff_t)tables.ball_end(det);
    auto stop = std::upper_bound(first, last, two_t, [](pm::cumulative_time_int bound, pm::weight_int w) {
        return bound < (pm::cumulative_time_int)w;
    });
    uint64_t walked = (uint64_t)(stop - first);
    return stop == last ? walked : walked + 1;
}

/// The union-find `find`, with path halving. Every link points a node at a **smaller** index, so
/// halving preserves that invariant and the root of a set is always its minimum member — which is
/// what makes the component enumeration sorted without a second sort (§0 determinism).
uint32_t find_root(std::vector<uint32_t>& parent, uint32_t node) {
    while (parent[node] != node) {
        parent[node] = parent[parent[node]];
        node = parent[node];
    }
    return node;
}

/// The component buffers' total capacity, the twin of `total_capacity` above and used the same way:
/// compared before and after an analysis to notice that it had to allocate.
uint64_t component_capacity(const BallComponents& c) {
    return c.component_of.capacity() + c.component_index.capacity() + c.roots.capacity() + c.sizes.capacity() +
           c.member_offsets.capacity() + c.members.capacity() + c.local_index.capacity() + c.adj_offsets.capacity() +
           c.adj_target.capacity() + c.adj_weight.capacity() + c.dijkstra_dist.capacity() + c.dijkstra_done.capacity() +
           c.bfs_dist.capacity() + c.bfs_queue.capacity() + c.fill_cursor.capacity();
}

/// The split buffers' total capacity, used exactly as the two above are: compared before and after
/// a decomposition to notice that it had to allocate.
uint64_t split_capacity(const BallComponentSplit& s) {
    return s.component_of.capacity() + s.component_index.capacity() + s.roots.capacity() + s.sizes.capacity() +
           s.member_offsets.capacity() + s.members.capacity() + s.local_index.capacity() + s.edge_offsets.capacity() +
           s.edge_index.capacity() + s.node_boundary.capacity() + s.status.capacity() + s.fill_cursor.capacity() +
           s.sub.h_to_det.capacity() + s.sub.edges.capacity() + s.sub.boundary_edges.capacity();
}

}  // namespace

void BallComponents::clear() {
    component_of.clear();
    component_index.clear();
    roots.clear();
    sizes.clear();
    member_offsets.clear();
    members.clear();
    local_index.clear();
    adj_offsets.clear();
    adj_target.clear();
    adj_weight.clear();
    fill_cursor.clear();
}

void BallComponentSplit::clear() {
    component_of.clear();
    component_index.clear();
    roots.clear();
    sizes.clear();
    member_offsets.clear();
    members.clear();
    local_index.clear();
    edge_offsets.clear();
    edge_index.clear();
    node_boundary.clear();
    status.clear();
    fill_cursor.clear();
    sub.clear();
}

void BallGraphArena::reset_for_graph(size_t num_detector_nodes) {
    graph.clear();
    det_to_h.assign(num_detector_nodes, NOT_A_DEFECT);
    syndrome_words.assign((num_detector_nodes + 63) / 64, 0);
    touched_words.clear();
    components.clear();
    components.grow_events = 0;
    split.clear();
    split.grow_events = 0;
    grow_events = 0;
}

void build_ball_graph(
    const BallTables& tables,
    const std::vector<uint64_t>& seeded_dets,
    horizon_int horizon,
    BallGraphArena& arena,
    BallGraphBuildMode mode,
    BallGraphTiming* timing,
    BallGraphCounts* counts) {
    if (horizon != pm::NO_HORIZON && horizon > tables.t_max_int)
        throw std::invalid_argument(
            "The requested horizon exceeds the compiled BallParams::T_max; recompile the ball tables. Decoding "
            "beyond T_max would silently drop reachable pairs.");
    if (arena.det_to_h.size() != tables.num_nodes)
        arena.reset_for_graph(tables.num_nodes);

    uint64_t capacity_before = total_capacity(arena);
    BallGraph& graph = arena.graph;
    graph.clear();

    // The horizon is finite in every truncated run; `NO_HORIZON` would ask for the complete defect
    // metric, which the ball tables do not hold and which `R >= 2 * T_max` was never sized for.
    if (horizon == pm::NO_HORIZON)
        throw std::invalid_argument(
            "The ball graph front end needs a finite horizon: at T = infinity `H` is the complete defect graph, "
            "which no ball radius covers. Decode that case on G with M1's path.");
    pm::cumulative_time_int two_t = 2 * (pm::cumulative_time_int)horizon;

    HiResTimer timer;
    if (timing != nullptr)
        timer.start();

    // ---- Nodes, and the two lookup structures the intersection reads.
    graph.h_to_det.reserve(seeded_dets.size());
    for (uint64_t det : seeded_dets) {
        if (det >= tables.num_nodes)
            throw std::invalid_argument("Detection event index is outside the detector graph.");
        assert(
            (graph.h_to_det.empty() || det > graph.h_to_det.back()) &&
            "seeded detection events must be sorted ascending and deduplicated");
        arena.det_to_h[det] = (uint32_t)graph.h_to_det.size();
        graph.h_to_det.push_back(det);

        size_t word = det / 64;
        if (arena.syndrome_words[word] == 0)
            arena.touched_words.push_back((uint32_t)word);
        arena.syndrome_words[word] |= (uint64_t)1 << (det % 64);
    }

    // ---- Edges. Each pair is emitted once, from its lower-id endpoint.
    //
    // The `structural_*` locals are the §M2 hardware-budget counters. Everything that is a plain
    // register increment is accumulated unconditionally and stored out once at the end, so the
    // no-profile path pays a loop-carried add and no memory traffic; anything that needs a table
    // read the decode itself does not do sits behind `counts != nullptr`.
    uint64_t structural_scan_entries = 0;
    uint64_t structural_scan_words = 0;
    uint64_t structural_hit_bytes = 0;
    uint64_t structural_other_mode_bytes = 0;
    uint64_t structural_edges = 0;
    uint64_t structural_boundary_edges = 0;
    for (uint32_t i = 0; i < graph.h_to_det.size(); i++) {
        uint64_t det = graph.h_to_det[i];

        if (mode == BallGraphBuildMode::SCAN) {
            uint64_t end = tables.ball_end(det);
            for (uint64_t e = tables.ball_begin(det); e < end; e++) {
                structural_scan_entries++;
                // Entries are sorted by `(w_int, target)`, so the first entry past `2T` ends the
                // scan: everything after it is unreachable by §M2.0, not merely unlikely.
                if ((pm::cumulative_time_int)tables.ball_w_int[e] > two_t)
                    break;
                uint32_t target = tables.ball_target[e];
                if (target <= det)
                    continue;
                uint32_t j = arena.det_to_h[target];
                if (j == NOT_A_DEFECT)
                    continue;
                graph.edges.push_back(BallGraphEdge{i, j, tables.ball_w_int[e], e});
                structural_edges++;
                if (counts != nullptr)
                    structural_hit_bytes += mask_bytes(tables, e);
            }
        } else {
            uint64_t words_begin = tables.ball_word_offsets[det];
            uint64_t words_end = tables.ball_word_offsets[det + 1];
            uint32_t base = tables.ball_word_base[det];
            uint64_t entry_begin = tables.ball_begin(det);
            for (uint64_t k = words_begin; k < words_end; k++) {
                size_t syndrome_word = base + (k - words_begin);
                if (syndrome_word >= arena.syndrome_words.size())
                    break;
                structural_scan_words++;
                uint64_t hits = tables.ball_words[k] & arena.syndrome_words[syndrome_word];
                while (hits != 0) {
                    uint64_t bit = hits & (~hits + 1);
                    int index = __builtin_ctzll(hits);
                    hits ^= bit;
                    uint32_t target = (uint32_t)(syndrome_word * 64 + index);
                    if (target <= det)
                        continue;
                    // Rank within the node's bitset is the index into the ascending-target entry
                    // pool, which is how the bitset path names the same canonical entry SCAN does.
                    uint32_t rank =
                        tables.ball_word_rank[k] + (uint32_t)__builtin_popcountll(tables.ball_words[k] & (bit - 1));
                    uint64_t entry = tables.ball_entry_by_rank[entry_begin + rank];
                    pm::weight_int w_int = tables.ball_w_int[entry];
                    // The rank indirection and the weight are read for every candidate, including
                    // the ones the `2T` test then rejects: the bitset carries neither, so there is
                    // no way to apply the test without fetching them first.
                    if (counts != nullptr)
                        structural_hit_bytes += ball_element_bytes::WORD_RANK + ball_element_bytes::ENTRY_BY_RANK +
                                                ball_element_bytes::WEIGHT;
                    if ((pm::cumulative_time_int)w_int > two_t)
                        continue;
                    graph.edges.push_back(BallGraphEdge{i, arena.det_to_h[target], w_int, entry});
                    structural_edges++;
                    if (counts != nullptr)
                        structural_hit_bytes += mask_bytes(tables, entry);
                }
            }
        }

        bool has_boundary = tables.has_bcost[det] != 0;
        bool boundary_within_horizon =
            has_boundary && (pm::cumulative_time_int)tables.bcost_w_int[det] <= (pm::cumulative_time_int)horizon;
        if (boundary_within_horizon) {
            graph.boundary_edges.push_back(BallBoundaryEdge{i, tables.bcost_w_int[det], det});
            structural_boundary_edges++;
        }
        if (counts != nullptr) {
            structural_hit_bytes += ball_element_bytes::HAS_BCOST;
            if (has_boundary)
                structural_hit_bytes += ball_element_bytes::BCOST_WEIGHT;
            if (boundary_within_horizon)
                structural_hit_bytes += boundary_mask_bytes(tables, det);
            structural_other_mode_bytes +=
                mode == BallGraphBuildMode::SCAN
                    ? bitset_words_counterfactual(tables, det, arena.syndrome_words.size()) * ball_element_bytes::WORD
                    : scan_entries_counterfactual(tables, det, two_t) *
                          (ball_element_bytes::TARGET + ball_element_bytes::WEIGHT);
        }
    }

    // The edge records emitted are exactly what `H` ends up holding: `graph.edges` is cleared at
    // the top and only ever appended to, and each undirected pair is appended once, from its
    // lower-id endpoint. That is what lets the aggregate read `hbld_edges_written` off `H`'s own
    // sizes and still call it a count of writes.
    assert(structural_edges == graph.edges.size() && "an undirected pair was emitted more than once");
    assert(structural_boundary_edges == graph.boundary_edges.size());

    if (counts != nullptr) {
        counts->isect_scan_bytes = structural_scan_entries * (ball_element_bytes::TARGET + ball_element_bytes::WEIGHT) +
                                   structural_scan_words * ball_element_bytes::WORD;
        counts->isect_hit_bytes = structural_hit_bytes;
        counts->isect_scan_bytes_other_mode = structural_other_mode_bytes;
        counts->edges_written = structural_edges;
        counts->boundary_edges_written = structural_boundary_edges;
    }

    if (timing != nullptr) {
        timing->intersect_ns = timer.elapsed_ns();
        timer.start();
    }

    // SCAN emits a node's hits in weight order and BITSET in target order. Sorting here is what
    // makes the two modes produce the *same* `H`, and what keeps `H` a function of the shot alone.
    std::sort(graph.edges.begin(), graph.edges.end(), [](const BallGraphEdge& a, const BallGraphEdge& b) {
        return a.i != b.i ? a.i < b.i : a.j < b.j;
    });

    // ---- Undo the syndrome bitset through the touched list; `det_to_h` through the node list.
    for (uint32_t word : arena.touched_words)
        arena.syndrome_words[word] = 0;
    arena.touched_words.clear();
    for (uint64_t det : graph.h_to_det)
        arena.det_to_h[det] = NOT_A_DEFECT;

    if (total_capacity(arena) != capacity_before)
        arena.grow_events++;
    if (timing != nullptr)
        timing->finalize_ns = timer.elapsed_ns();
}

void decompose_ball_components(const BallGraph& graph, BallComponentSplit& s) {
    uint64_t capacity_before = split_capacity(s);
    uint32_t n = (uint32_t)graph.num_nodes();

    // ---- Union-find over the defect-defect edges, and only those. A boundary edge joins nothing:
    // the boundary is not a node of `H`, so two defects that both reach it are not thereby
    // connected, and treating them as connected would merge components that never interact.
    s.component_of.resize(n);
    for (uint32_t i = 0; i < n; i++)
        s.component_of[i] = i;
    for (const BallGraphEdge& edge : graph.edges) {
        uint32_t a = find_root(s.component_of, edge.i);
        uint32_t b = find_root(s.component_of, edge.j);
        if (a == b)
            continue;
        // Always link the larger root under the smaller, so a set's root is its minimum member.
        if (a < b)
            s.component_of[b] = a;
        else
            s.component_of[a] = b;
    }
    // Every non-root points at a strictly smaller index, so one ascending pass leaves every entry
    // pointing straight at its root — no second `find` walk, and no recursion.
    for (uint32_t i = 0; i < n; i++) {
        uint32_t parent = s.component_of[i];
        s.component_of[i] = parent == i ? i : s.component_of[parent];
    }

    // ---- Components, in ascending root order, and their sizes.
    s.component_index.resize(n);
    s.roots.clear();
    s.sizes.clear();
    for (uint32_t i = 0; i < n; i++) {
        if (s.component_of[i] == i) {
            s.component_index[i] = (uint32_t)s.roots.size();
            s.roots.push_back(i);
            s.sizes.push_back(0);
        }
    }
    for (uint32_t i = 0; i < n; i++)
        s.sizes[s.component_index[s.component_of[i]]]++;

    // ---- Members, grouped by component and ascending within each. `i` ascends, so the block of a
    // component is filled in ascending `H`-node order for free — which is §1's bit-exactness
    // clause: the sub-`H` relabelling preserves the order the nodes have in the full `H`, so
    // blossom's equal-time tie-break inside the component is unchanged.
    size_t num_components = s.roots.size();
    s.member_offsets.assign(num_components + 1, 0);
    for (size_t c = 0; c < num_components; c++)
        s.member_offsets[c + 1] = s.member_offsets[c] + s.sizes[c];
    s.members.resize(n);
    s.local_index.resize(n);
    s.fill_cursor.assign(s.member_offsets.begin(), s.member_offsets.end());
    for (uint32_t i = 0; i < n; i++) {
        uint32_t index = s.component_index[s.component_of[i]];
        uint32_t slot = s.fill_cursor[index]++;
        s.members[slot] = i;
        s.local_index[i] = slot - s.member_offsets[index];
    }

    // ---- Internal edges, grouped the same way. Both endpoints of an `H` edge share a component,
    // so every edge lands in exactly one block, and `e` ascends so each block is ascending in the
    // global edge order — which `graph.edges` holds sorted by `(i, j)`.
    s.edge_offsets.assign(num_components + 1, 0);
    for (const BallGraphEdge& edge : graph.edges)
        s.edge_offsets[s.component_index[s.component_of[edge.i]] + 1]++;
    for (size_t c = 0; c < num_components; c++)
        s.edge_offsets[c + 1] += s.edge_offsets[c];
    s.edge_index.resize(graph.edges.size());
    s.fill_cursor.assign(s.edge_offsets.begin(), s.edge_offsets.end());
    for (uint32_t e = 0; e < (uint32_t)graph.edges.size(); e++) {
        const BallGraphEdge& edge = graph.edges[e];
        assert(
            s.component_of[edge.i] == s.component_of[edge.j] &&
            "an H edge joins two nodes the union-find left in different components");
        s.edge_index[s.fill_cursor[s.component_index[s.component_of[edge.i]]]++] = e;
    }

    // ---- Where each node's boundary edge is, if it has one. `H` gives a node at most one, so this
    // is a single slot; §2.4 copies `H`'s own boundary edges rather than re-testing `bcost <= T`.
    s.node_boundary.assign(n, BallComponentSplit::NO_BOUNDARY);
    for (uint32_t b = 0; b < (uint32_t)graph.boundary_edges.size(); b++) {
        uint32_t node = graph.boundary_edges[b].i;
        assert(s.node_boundary[node] == BallComponentSplit::NO_BOUNDARY && "H gave a node two boundary edges");
        s.node_boundary[node] = b;
    }

    s.status.assign(num_components, BallComponentSplit::COMPONENT_COMPLETE);

    if (split_capacity(s) != capacity_before)
        s.grow_events++;
}

void build_component_subgraph(const BallGraph& graph, BallComponentSplit& s, size_t index) {
    uint64_t capacity_before = split_capacity(s);
    BallGraph& sub = s.sub;
    sub.clear();

    uint32_t begin = s.member_offsets[index];
    uint32_t size = s.sizes[index];

    // Local node `j` <-> `members[begin + j]`, and the members are ascending, so `sub.h_to_det` is a
    // strictly ascending subsequence of `graph.h_to_det`: a residual sorted here is sorted in `G`,
    // exactly as it is when the solve runs on the whole of `H`.
    sub.h_to_det.reserve(size);
    for (uint32_t j = 0; j < size; j++)
        sub.h_to_det.push_back(graph.h_to_det[s.members[begin + j]]);

    // Weights and observable entries are **copied** from `H`, so extraction is bit-identical: the
    // same integers and the same canonical ball entries the monolithic instance would have used.
    // The block is ascending in the global edge order and `local_index` is monotone within a
    // component, so `sub.edges` inherits `H`'s `(i, j)` sort with `i < j` — which
    // `BallMwpm::rebuild` relies on for its ascending-neighbour adjacency.
    for (uint32_t k = s.edge_offsets[index]; k < s.edge_offsets[index + 1]; k++) {
        const BallGraphEdge& edge = graph.edges[s.edge_index[k]];
        sub.edges.push_back(
            BallGraphEdge{s.local_index[edge.i], s.local_index[edge.j], edge.w_int, edge.entry});
    }

    // Exactly the boundary edges `H` already has for these nodes — the cutoff is not re-derived
    // (§2.4). Members ascend, so the list comes out sorted by local index, which is what
    // `BallMwpm::rebuild` requires of `boundary_edges`.
    for (uint32_t j = 0; j < size; j++) {
        uint32_t b = s.node_boundary[s.members[begin + j]];
        if (b == BallComponentSplit::NO_BOUNDARY)
            continue;
        const BallBoundaryEdge& edge = graph.boundary_edges[b];
        sub.boundary_edges.push_back(BallBoundaryEdge{j, edge.w_int, edge.det});
    }

    if (split_capacity(s) != capacity_before)
        s.grow_events++;
}

void analyze_ball_components(const BallGraph& graph, BallComponents& c, uint32_t diameter_cap) {
    uint64_t capacity_before = component_capacity(c);
    uint32_t n = (uint32_t)graph.num_nodes();

    // ---- Union-find over the defect-defect edges, and only those. A boundary edge joins nothing:
    // the boundary is not a node of `H`, so two defects that both reach it are not thereby
    // connected, and treating them as connected would merge components the solver keeps apart.
    c.component_of.resize(n);
    for (uint32_t i = 0; i < n; i++)
        c.component_of[i] = i;
    for (const BallGraphEdge& edge : graph.edges) {
        uint32_t a = find_root(c.component_of, edge.i);
        uint32_t b = find_root(c.component_of, edge.j);
        if (a == b)
            continue;
        // Always link the larger root under the smaller, so a set's root is its minimum member.
        if (a < b)
            c.component_of[b] = a;
        else
            c.component_of[a] = b;
    }
    // Every non-root points at a strictly smaller index, so one ascending pass leaves every entry
    // pointing straight at its root — no second `find` walk, and no recursion.
    for (uint32_t i = 0; i < n; i++) {
        uint32_t parent = c.component_of[i];
        c.component_of[i] = parent == i ? i : c.component_of[parent];
    }

    // ---- Components, in ascending root order, and their sizes.
    c.component_index.resize(n);
    c.roots.clear();
    c.sizes.clear();
    for (uint32_t i = 0; i < n; i++) {
        if (c.component_of[i] == i) {
            c.component_index[i] = (uint32_t)c.roots.size();
            c.roots.push_back(i);
            c.sizes.push_back(0);
        }
    }
    for (uint32_t i = 0; i < n; i++)
        c.sizes[c.component_index[c.component_of[i]]]++;

    // ---- Members, grouped by component and ascending within each. `i` ascends, so the block of a
    // component is filled in ascending `H`-node order for free.
    size_t num_components = c.roots.size();
    c.member_offsets.assign(num_components + 1, 0);
    for (size_t k = 0; k < num_components; k++)
        c.member_offsets[k + 1] = c.member_offsets[k] + c.sizes[k];
    c.members.resize(n);
    c.local_index.resize(n);
    c.fill_cursor.assign(c.member_offsets.begin(), c.member_offsets.end());
    for (uint32_t i = 0; i < n; i++) {
        uint32_t index = c.component_index[c.component_of[i]];
        uint32_t slot = c.fill_cursor[index]++;
        c.members[slot] = i;
        c.local_index[i] = slot - c.member_offsets[index];
    }

    // ---- Undirected adjacency, both directions of every defect-defect edge.
    c.adj_offsets.assign((size_t)n + 1, 0);
    for (const BallGraphEdge& edge : graph.edges) {
        c.adj_offsets[edge.i + 1]++;
        c.adj_offsets[edge.j + 1]++;
    }
    for (uint32_t i = 0; i < n; i++)
        c.adj_offsets[i + 1] += c.adj_offsets[i];
    c.adj_target.resize(2 * graph.edges.size());
    c.adj_weight.resize(2 * graph.edges.size());
    c.fill_cursor.assign(c.adj_offsets.begin(), c.adj_offsets.end());
    for (const BallGraphEdge& edge : graph.edges) {
        uint32_t slot = c.fill_cursor[edge.i]++;
        c.adj_target[slot] = edge.j;
        c.adj_weight[slot] = edge.w_int;
        slot = c.fill_cursor[edge.j]++;
        c.adj_target[slot] = edge.i;
        c.adj_weight[slot] = edge.w_int;
    }

    c.dijkstra_dist.resize(diameter_cap);
    c.dijkstra_done.resize(diameter_cap);
    c.bfs_dist.resize(diameter_cap);
    c.bfs_queue.resize(diameter_cap);

    if (component_capacity(c) != capacity_before)
        c.grow_events++;
}

pm::cumulative_time_int component_diameter(BallComponents& c, size_t index, uint32_t diameter_cap) {
    uint32_t size = c.sizes[index];
    if (size > diameter_cap || size > c.dijkstra_dist.size())
        return -1;
    if (size <= 1)
        return 0;

    constexpr pm::cumulative_time_int UNREACHED = std::numeric_limits<pm::cumulative_time_int>::max();
    uint32_t begin = c.member_offsets[index];
    pm::cumulative_time_int diameter = 0;
    for (uint32_t source = 0; source < size; source++) {
        for (uint32_t k = 0; k < size; k++) {
            c.dijkstra_dist[k] = UNREACHED;
            c.dijkstra_done[k] = 0;
        }
        c.dijkstra_dist[source] = 0;
        for (uint32_t settled = 0; settled < size; settled++) {
            // A linear scan for the minimum rather than a heap: the component is capped at
            // `MAX_DIAMETER_COMPONENT_SIZE` members, where the scan is faster and, unlike a
            // priority queue, allocates nothing.
            uint32_t best = size;
            pm::cumulative_time_int best_dist = UNREACHED;
            for (uint32_t k = 0; k < size; k++) {
                if (!c.dijkstra_done[k] && c.dijkstra_dist[k] < best_dist) {
                    best_dist = c.dijkstra_dist[k];
                    best = k;
                }
            }
            // Cannot happen — the members of a component are connected through it by definition —
            // but stopping on it is cheaper than trusting it.
            if (best == size)
                break;
            c.dijkstra_done[best] = 1;
            uint32_t node = c.members[begin + best];
            for (uint32_t e = c.adj_offsets[node]; e < c.adj_offsets[node + 1]; e++) {
                uint32_t neighbour = c.local_index[c.adj_target[e]];
                pm::cumulative_time_int candidate = best_dist + (pm::cumulative_time_int)c.adj_weight[e];
                if (candidate < c.dijkstra_dist[neighbour])
                    c.dijkstra_dist[neighbour] = candidate;
            }
        }
        for (uint32_t k = 0; k < size; k++) {
            if (c.dijkstra_dist[k] != UNREACHED)
                diameter = std::max(diameter, c.dijkstra_dist[k]);
        }
    }
    return diameter;
}

int32_t component_hop_diameter(BallComponents& c, size_t index, uint32_t diameter_cap) {
    uint32_t size = c.sizes[index];
    if (size > diameter_cap || size > c.bfs_dist.size())
        return -1;
    if (size <= 1)
        return 0;

    uint32_t begin = c.member_offsets[index];
    int32_t diameter = 0;
    for (uint32_t source = 0; source < size; source++) {
        for (uint32_t k = 0; k < size; k++)
            c.bfs_dist[k] = -1;
        c.bfs_dist[source] = 0;
        c.bfs_queue[0] = source;
        uint32_t head = 0;
        uint32_t tail = 1;
        while (head < tail) {
            uint32_t local = c.bfs_queue[head++];
            uint32_t node = c.members[begin + local];
            for (uint32_t e = c.adj_offsets[node]; e < c.adj_offsets[node + 1]; e++) {
                uint32_t neighbour = c.local_index[c.adj_target[e]];
                if (c.bfs_dist[neighbour] >= 0)
                    continue;
                c.bfs_dist[neighbour] = c.bfs_dist[local] + 1;
                diameter = std::max(diameter, c.bfs_dist[neighbour]);
                c.bfs_queue[tail++] = neighbour;
            }
        }
    }
    return diameter;
}

}  // namespace spec_matching
}  // namespace pm
