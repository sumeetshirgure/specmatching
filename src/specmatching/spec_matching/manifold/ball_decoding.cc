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

#include "specmatching/spec_matching/manifold/ball_decoding.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <sstream>

#include "specmatching/spec_matching/certificate/max_dual.h"
#include "specmatching/spec_matching/manifold/ball_serialize.h"

namespace pm {
namespace spec_matching {

namespace {

std::string describe(const std::vector<uint64_t>& values) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < values.size(); i++)
        out << (i ? ", " : "") << values[i];
    out << "]";
    return out.str();
}

/// `bcost_int(x)`, read from the ball tables and nowhere else. `exists == false` is the design's
/// `+inf`: the defect has no boundary path within `R`, so `H` gave it no boundary edge and no
/// boundary match of it is legal.
struct BoundaryCost {
    bool exists{false};
    pm::cumulative_time_int w_int{0};
};

BoundaryCost boundary_cost(const BallTables& tables, uint64_t det) {
    if (tables.has_bcost[det] == 0)
        return BoundaryCost();
    return BoundaryCost{true, (pm::cumulative_time_int)tables.bcost_w_int[det]};
}

/// What the small-component resolver did with a component. `COMMIT` means it settled the component
/// exactly and performed the matches; `RESIDUAL` means the component has no feasible matching at
/// all and the shot escalates on it; `SOLVER` means it was larger than `k` and was handed on.
///
/// There is no fourth outcome. Ties are broken here, not deferred, and an exact-`T` match is legal
/// rather than ambiguous — nothing routes to the solver except by being size `> k`.
enum class SmallVerdict : uint8_t { SOLVER = 0, COMMIT = 1, RESIDUAL = 2 };

constexpr uint32_t MAX_SMALL = BallPrune::MAX_SMALL_COMPONENT_SIZE;

/// The resolver's view of one component of size 1..4: which matches are legal, and what they weigh,
/// in integer time units and in local indices `0..size-1` (ascending `H`-node order).
///
/// Legality comes from **`H`'s own inclusion cutoffs**, so the resolver and `H` agree by
/// construction: a defect-defect pairing is legal iff `H` gave the pair an edge, i.e.
/// `d_G(u,v) <= 2T`, and a boundary match is legal iff `bcost_int(x) <= T_int`. Both are weak
/// inequalities, exactly as `H`'s are; there is no separate exact-`T` event to route around.
struct SmallProblem {
    uint32_t size{0};
    bool boundary_legal[MAX_SMALL]{};
    pm::cumulative_time_int boundary_w[MAX_SMALL]{};
    bool pair_legal[MAX_SMALL][MAX_SMALL]{};
    pm::cumulative_time_int pair_w[MAX_SMALL][MAX_SMALL]{};
};

/// One configuration of "defect-pairings-with-boundary-fill": a set of disjoint legal
/// defect-defect pairs, with every unpaired member taking a legal boundary match.
struct SmallMatching {
    /// False on the empty result, i.e. the component admits no such configuration at all.
    bool feasible{false};
    pm::cumulative_time_int weight{0};
    uint32_t num_pairs{0};
    /// Local indices of the paired members, `pair_a[i] < pair_b[i]`, and the list ascending in
    /// `pair_a` — which is the sorted `(min, max)` pair list the tie-break compares.
    uint32_t pair_a[MAX_SMALL / 2]{};
    uint32_t pair_b[MAX_SMALL / 2]{};
    uint32_t num_boundary{0};
    uint32_t boundary[MAX_SMALL]{};
};

/// The order the resolver minimises in: total `w_int` first, and among equal-weight configurations
/// the lexicographically smallest sorted `(min, max)` defect pair list, with a shorter list — one
/// that is a prefix of the other — ordered first.
///
/// Equal weight means equal committed weight, and an observable that differs by at most a
/// homologically trivial cycle (zero flip), so either configuration is exact. The tie-break decides
/// only *which* one, so that two runs are bit-identical (§0).
bool small_matching_better(const SmallMatching& a, const SmallMatching& b) {
    if (!b.feasible)
        return a.feasible;
    if (!a.feasible)
        return false;
    if (a.weight != b.weight)
        return a.weight < b.weight;
    for (uint32_t i = 0; i < a.num_pairs && i < b.num_pairs; i++) {
        if (a.pair_a[i] != b.pair_a[i])
            return a.pair_a[i] < b.pair_a[i];
        if (a.pair_b[i] != b.pair_b[i])
            return a.pair_b[i] < b.pair_b[i];
    }
    return a.num_pairs < b.num_pairs;
}

/// Extends `current` from the lowest unmatched member and keeps the best complete configuration.
///
/// Taking the *lowest* unmatched member at every step enumerates each configuration exactly once,
/// and emits its pairs in ascending `pair_a` order, which is the sorted pair list the tie-break
/// wants without a sort. The candidate sets of the design's sizes 1–4 are what this generates:
/// a singleton's boundary match; a pair's `u–v` or `u,v` both to boundary; a triple's three
/// "one to boundary, the other two paired" plus all-three-to-boundary; a quadruple's three internal
/// perfect matchings, its `C(4,2)` "one pair plus two boundaries", and all-four-to-boundary.
void extend_small_matching(const SmallProblem& problem, uint32_t used, SmallMatching& current, SmallMatching& best) {
    uint32_t i = 0;
    while (i < problem.size && (used & (1u << i)) != 0)
        i++;
    if (i == problem.size) {
        current.feasible = true;
        if (small_matching_better(current, best))
            best = current;
        return;
    }

    // The boundary fill is tried before the pairings, so that among equal-weight configurations the
    // shorter pair list is reached first. That is the tie-break's order made explicit in the
    // enumeration; `small_matching_better` still decides, so the two cannot drift apart.
    if (problem.boundary_legal[i]) {
        current.boundary[current.num_boundary++] = i;
        current.weight += problem.boundary_w[i];
        extend_small_matching(problem, used | (1u << i), current, best);
        current.weight -= problem.boundary_w[i];
        current.num_boundary--;
    }
    for (uint32_t j = i + 1; j < problem.size; j++) {
        if ((used & (1u << j)) != 0 || !problem.pair_legal[i][j])
            continue;
        current.pair_a[current.num_pairs] = i;
        current.pair_b[current.num_pairs] = j;
        current.num_pairs++;
        current.weight += problem.pair_w[i][j];
        extend_small_matching(problem, used | (1u << i) | (1u << j), current, best);
        current.weight -= problem.pair_w[i][j];
        current.num_pairs--;
    }
}

/// Sizes 1–4, exactly: the minimum-weight pairing-with-boundary-fill, or infeasible.
///
/// Infeasible is not a failure and not an error. Even size does not guarantee a complete matching —
/// a size-4 "star" whose centre is within `2T` of three mutually far members, none of them with a
/// legal boundary, has no feasible pairing, and neither does an odd component with no legal
/// boundary. The caller routes those to the residual and the shot escalates, exactly as the solver
/// would have left it.
SmallMatching resolve_small_component(const SmallProblem& problem) {
    SmallMatching best;
    SmallMatching current;
    extend_small_matching(problem, 0, current, best);
    return best;
}

/// Position of an `H`-node within its component's ascending member block, i.e. the local index the
/// two structs above are written in. A linear scan over at most four entries.
template <typename Members>
uint32_t local_index_of(const Members& members, uint32_t size, uint32_t node) {
    for (uint32_t a = 0; a < size; a++) {
        if (members[a] == node)
            return a;
    }
    assert(false && "an H edge names a node outside the component it is internal to");
    return 0;
}

/// The observable mask of one ball entry's stored path, XORed out of the entry's observable id
/// list. Zero above 64 observables, exactly as `BallMwpm::rebuild` leaves it there: the mask is
/// unusable and the match-edge flavour is the only correct readout (§0's observables rule).
pm::obs_int path_obs_mask(const BallTables& tables, uint64_t entry, bool use_masks) {
    if (!use_masks)
        return 0;
    pm::obs_int mask = 0;
    for (uint64_t k = tables.ball_mask_offsets[entry]; k < tables.ball_mask_offsets[entry + 1]; k++)
        mask ^= (pm::obs_int)1 << tables.ball_mask_ids[k];
    return mask;
}

pm::obs_int boundary_obs_mask(const BallTables& tables, uint64_t det, bool use_masks) {
    if (!use_masks)
        return 0;
    pm::obs_int mask = 0;
    for (uint64_t k = tables.bcost_mask_offsets[det]; k < tables.bcost_mask_offsets[det + 1]; k++)
        mask ^= (pm::obs_int)1 << tables.bcost_mask_ids[k];
    return mask;
}

void sort_pairs(std::vector<CommittedPair>& pairs) {
    for (auto& pair : pairs) {
        if (pair.to >= 0 && pair.to < pair.from)
            std::swap(pair.from, pair.to);
    }
    std::sort(pairs.begin(), pairs.end(), [](const CommittedPair& a, const CommittedPair& b) {
        return a.from != b.from ? a.from < b.from : a.to < b.to;
    });
}

std::string describe(const std::vector<CommittedPair>& pairs) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < pairs.size(); i++)
        out << (i ? ", " : "") << "(" << pairs[i].from << ", " << pairs[i].to << ")";
    out << "]";
    return out.str();
}

}  // namespace

BallDecoder BallDecoder::from_detector_error_model(
    const stim::DetectorErrorModel& dem,
    BallConfig config,
    pm::weight_int num_distinct_weights,
    const char* ball_artifact_path) {
    // The search flooder is only needed when the obs_int masks are unusable, which is exactly when
    // the ball tables have to read observable ids off the search graph instead.
    bool need_search_graph = dem.count_observables() > sizeof(pm::obs_int) * 8;
    pm::Mwpm mwpm = pm::detector_error_model_to_mwpm(dem, num_distinct_weights, need_search_graph);
    return from_mwpm(std::move(mwpm), std::move(config), ball_artifact_path);
}

BallDecoder BallDecoder::from_mwpm(pm::Mwpm g_mwpm, BallConfig config, const char* ball_artifact_path) {
    BallDecoder decoder;
    decoder.config = std::move(config);
    decoder.g_mwpm = std::move(g_mwpm);
    decoder.finish_construction(ball_artifact_path);
    return decoder;
}

void BallDecoder::finish_construction(const char* ball_artifact_path) {
    config.ball.validate();
    if (config.T > config.ball.T_max)
        throw std::invalid_argument(
            "BallConfig::T exceeds BallParams::T_max; the compiled tables do not cover that horizon.");
    // §M7: the M1 oracle does not apply on the stock path, because that path produces no truncated
    // intermediate state to reproduce. Its oracle is stock exact decode on `G` (§M7.6 level 1),
    // which the M7 test suite runs; silently accepting the flag and comparing against a harvest
    // that never happens would be worse than refusing it.
    if (config.stock_on_h && config.verify_against_g)
        throw std::invalid_argument(
            "verify_against_g compares H's truncated harvest against M1 on G; the stock-on-H path (§M7) has no "
            "truncated harvest. Its oracle is stock exact decode on G — see §M7.6 level 1.");
    // §A. `k` is capped at 4 because no resolver is defined above it: the enumeration of
    // pairings-with-boundary-fill is written out for sizes 1–4 and for nothing larger. Rejected here
    // rather than silently clamped, so a run cannot report a `k` it did not decode at.
    if (config.prune_component_max_size < 0 ||
        config.prune_component_max_size > (int)BallPrune::MAX_SMALL_COMPONENT_SIZE)
        throw std::invalid_argument(
            "prune_component_max_size (k) must be in [0, " + std::to_string(BallPrune::MAX_SMALL_COMPONENT_SIZE) +
            "]; no small-component resolver is defined above that.");
    assert(config.prune_component_max_size <= (int)BallPrune::MAX_SMALL_COMPONENT_SIZE && "k exceeds the resolver cap");

    const pm::MatchingGraph& graph = g_mwpm.flooder.graph;
    horizon = to_time_units(config.T, graph.normalising_constant);

    if (ball_artifact_path != nullptr) {
        tables = load_ball_tables(ball_artifact_path, graph, config.ball);
    } else {
        tables = compile_ball_tables(g_mwpm, config.ball, config.compile_threads);
    }
    // Invariant 6, restated at decoder construction: `R >= 2 * T_max` in the units that actually
    // matter, and one shared normalising constant between `G` and `H` (§0 unit rule).
    if (tables.normalising_constant != graph.normalising_constant)
        throw std::invalid_argument("Ball tables were compiled with a different normalising constant.");
    if (tables.r_int < 2 * tables.t_max_int)
        throw std::invalid_argument("Ball tables violate R >= 2 * T_max in time units.");
    if (horizon > tables.t_max_int)
        throw std::invalid_argument("BallConfig::T exceeds the compiled T_max in time units.");

    // §B. One look at `G`'s own record of its negative-weight edges, taken once here rather than
    // per shot. The discretised weight sum is non-zero iff some edge came in negative; the event
    // set is what the §M2.1 preamble would symmetric-difference in, and it can cancel to empty
    // while the sum does not (the set is toggled per endpoint), so both are read.
    dem_has_negative_weights =
        g_mwpm.flooder.negative_weight_sum != 0 || !g_mwpm.flooder.negative_weight_detection_events.empty();

    h_mwpm.configure(graph.num_observables, graph.normalising_constant);
    arena.reset_for_graph(graph.nodes.size());
}

const BallGraph& BallDecoder::resolve_small_components(
    const BallGraph& h, pm::MatchingResult& resolved, SmallCommitCounts& counts) {
    // Untimed by construction: no timer is started here and none may be added (hard constraint 1).
    // This is a serial pre-pass on the critical path, excluded from what this branch reports by
    // measurement scope; see `decode_impl`'s call site.
    BallPrune& prune = arena.prune;
    compute_prune_components(h, prune);

    // The horizon in the decoder's own stored integer units. Read, never re-derived: converting
    // `config.T` a second time here would be the §0 unit-rule trap.
    const pm::cumulative_time_int t_int = horizon;
    bool use_masks = g_mwpm.flooder.graph.num_observables <= sizeof(pm::obs_int) * 8;
    // `k`, capped at the largest size a resolver is defined for. The cap is enforced at
    // construction; clamping again here is what keeps this loop's fixed-width blocks in range no
    // matter how the config was reached.
    uint32_t max_size = config.prune_component_max_size <= 0
                            ? 0
                            : std::min<uint32_t>((uint32_t)config.prune_component_max_size, MAX_SMALL);
    uint32_t n = (uint32_t)h.num_nodes();

    auto commit_boundary = [&](uint32_t node) {
        uint64_t det = h.h_to_det[node];
        resolved.obs_mask ^= boundary_obs_mask(tables, det, use_masks);
        resolved.weight += (pm::total_weight_int)tables.bcost_w_int[det];
        resolved_pairs.push_back(CommittedPair{(int64_t)det, -1, CommittedPair::NO_BALL_ENTRY});
        counts.boundary++;
    };

    // ---- One pass over the components, in ascending root order (§0 determinism). Roots are the
    // minimum member of their set, so `root` ascending is `H`-node order and the accumulation below
    // is a function of `H` alone.
    for (uint32_t root = 0; root < n; root++) {
        if (prune.component_of[root] != root)
            continue;
        uint32_t size = prune.component_size[root];
        if (size > max_size) {
            prune.verdict[root] = (uint8_t)SmallVerdict::SOLVER;
            continue;
        }

        // The component's legality table, off `H` and the ball tables and nothing else. The
        // presence of an `H` edge *is* the pairing's legality test (`d_G <= 2T`), and its `w_int`
        // is the pairing's weight, so the resolver cannot disagree with `H` about either.
        const auto& members = prune.small_members[root];
        assert(prune.small_member_count[root] == size && "the small-component member block is short");
        SmallProblem problem;
        problem.size = size;
        for (uint32_t a = 0; a < size; a++) {
            BoundaryCost cost = boundary_cost(tables, h.h_to_det[members[a]]);
            problem.boundary_legal[a] = cost.exists && cost.w_int <= t_int;
            problem.boundary_w[a] = cost.w_int;
        }
        uint32_t edge_of[MAX_SMALL][MAX_SMALL] = {};
        for (uint32_t e = 0; e < prune.small_edge_count[root]; e++) {
            uint32_t index = prune.small_edges[root][e];
            const BallGraphEdge& edge = h.edges[index];
            uint32_t a = local_index_of(members, size, edge.i);
            uint32_t b = local_index_of(members, size, edge.j);
            problem.pair_legal[a][b] = problem.pair_legal[b][a] = true;
            problem.pair_w[a][b] = problem.pair_w[b][a] = (pm::cumulative_time_int)edge.w_int;
            edge_of[a][b] = edge_of[b][a] = index;
        }

        SmallMatching matching = resolve_small_component(problem);
        if (!matching.feasible) {
            // No configuration at all: the component survives past the horizon and the shot
            // escalates on it, exactly as the solver would have left it.
            counts.residual++;
            prune.verdict[root] = (uint8_t)SmallVerdict::RESIDUAL;
            continue;
        }
        for (uint32_t p = 0; p < matching.num_pairs; p++) {
            // Mask and weight come off the ball entry behind the edge, which is the same entry
            // harvest would have extracted for this match — that is what keeps the committed
            // observables and weight bit-exact.
            const BallGraphEdge& edge = h.edges[edge_of[matching.pair_a[p]][matching.pair_b[p]]];
            resolved.obs_mask ^= path_obs_mask(tables, edge.entry, use_masks);
            resolved.weight += (pm::total_weight_int)edge.w_int;
            resolved_pairs.push_back(
                CommittedPair{(int64_t)h.h_to_det[edge.i], (int64_t)h.h_to_det[edge.j], edge.entry});
            counts.pairs++;
        }
        for (uint32_t b = 0; b < matching.num_boundary; b++)
            commit_boundary(members[matching.boundary[b]]);
        counts.defects_resolved += (int)size;
        prune.verdict[root] = (uint8_t)SmallVerdict::COMMIT;
    }

    // ---- Spread each component's verdict to its members, and count what is left for the solver.
    uint32_t solver_nodes = 0;
    for (uint32_t i = 0; i < n; i++) {
        // A RESIDUAL component is neither committed nor solved: the shot escalates on it and §A.5
        // discards Phase 1 in full, so there is nothing for the solver to learn from it.
        bool to_solver = (SmallVerdict)prune.verdict[prune.component_of[i]] == SmallVerdict::SOLVER;
        prune.node_to_solver[i] = to_solver ? 1 : 0;
        solver_nodes += to_solver ? 1 : 0;
    }

    // Nothing was resolved away: hand the solver `H` itself rather than copying it into an
    // identical sub-graph. This is the shot class the prune cannot help, and it should not pay for
    // being looked at.
    if (solver_nodes == n) {
        assert(counts.pairs == 0 && counts.boundary == 0 && counts.residual == 0 && counts.defects_resolved == 0);
        return h;
    }
    build_solver_subgraph(h, prune);
    return prune.solver_graph;
}

void BallDecoder::save_ball_artifact(const std::string& path) const {
    save_ball_tables(tables, path);
}

void BallDecoder::analyze_last_shot_components(BallProfile& prof, ComponentHistograms& histograms) {
    // No timer is started anywhere in this function, and none may be added: §C is untimed by
    // construction and the caller has already closed the shot's window.
    ComponentStats& stats = prof.components;
    stats = ComponentStats();
    histograms.clear();

    const BallGraph& h = arena.graph;
    BallComponents& components = arena.components;
    analyze_ball_components(h, components);

    // The horizon in the decoder's own stored integer units. Read, never re-derived: converting
    // `config.T` a second time here would be the §0 unit-rule trap.
    const pm::cumulative_time_int t_int = horizon;
    uint32_t num_nodes = (uint32_t)h.num_nodes();
    // The same `k` the decode path resolves at, clamped the same way, so the classification below
    // describes the run it is collected beside rather than a fixed size-2 rule.
    uint32_t max_size = config.prune_component_max_size <= 0
                            ? 0
                            : std::min<uint32_t>((uint32_t)config.prune_component_max_size, MAX_SMALL);

    stats.measured = 1;
    stats.num_components = (int)components.num_components();

    // ---- Per-edge and per-node distributions. Every `H` edge weight is an exact `d_G` between two
    // defects, which is the path-length distribution §C.2 asks for.
    for (const BallGraphEdge& edge : h.edges) {
        histograms.edge_weight_hist[ComponentHistograms::weight_bin(
            (pm::cumulative_time_int)edge.w_int, t_int, ComponentHistograms::WEIGHT_HIST_BINS)]++;
    }
    for (uint32_t i = 0; i < num_nodes; i++) {
        histograms.degree_hist[ComponentHistograms::count_bin(
            components.degree_of(i), ComponentHistograms::DEGREE_HIST_BINS)]++;
        BoundaryCost cost = boundary_cost(tables, h.h_to_det[i]);
        if (cost.exists) {
            histograms.bcost_hist[ComponentHistograms::weight_bin(
                cost.w_int, t_int, ComponentHistograms::BCOST_HIST_BINS)]++;
        }
    }

    // ---- Per component, in ascending root order (§0 determinism).
    for (size_t c = 0; c < components.num_components(); c++) {
        uint32_t size = components.sizes[c];
        uint32_t begin = components.member_offsets[c];
        histograms.size_hist[ComponentHistograms::count_bin(size, ComponentHistograms::SIZE_HIST_BINS)]++;
        stats.largest_component_size = std::max(stats.largest_component_size, (int)size);

        pm::cumulative_time_int diameter = component_diameter(components, c);
        if (diameter < 0) {
            stats.diameter_uncomputed_components++;
        } else {
            histograms.diameter_hist[ComponentHistograms::weight_bin(
                diameter, t_int, ComponentHistograms::DIAMETER_HIST_BINS)]++;
            stats.max_component_diameter_wint = std::max(stats.max_component_diameter_wint, (int)diameter);
        }

        bool touches_boundary = false;
        for (uint32_t k = 0; k < size; k++) {
            BoundaryCost cost = boundary_cost(tables, h.h_to_det[components.members[begin + k]]);
            if (cost.exists && cost.w_int <= t_int) {
                touches_boundary = true;
                break;
            }
        }
        stats.num_boundary_touching_components += touches_boundary ? 1 : 0;

        // The size classes, which are structural and say nothing about `k`: "trivial" is size <= 2
        // in this table whatever the resolver is configured to attempt, so the §C series stays
        // comparable across runs.
        if (size == 1) {
            stats.num_singleton_components++;
        } else if (size == 2) {
            stats.num_pair_components++;
        } else {
            stats.num_nontrivial_components++;
        }
        if (size <= 2) {
            stats.num_trivial_components++;
            stats.defects_in_trivial_components += (int)size;
        }

        // The resolver's classification, at the `k` this decoder is configured with. Nothing is
        // committed here: the verdict is counted and discarded, so this is a measurement of what
        // the prune removes and not the prune itself.
        SmallVerdict verdict = SmallVerdict::SOLVER;
        if (size <= max_size) {
            SmallProblem problem;
            problem.size = size;
            for (uint32_t a = 0; a < size; a++) {
                uint32_t node = components.members[begin + a];
                BoundaryCost cost = boundary_cost(tables, h.h_to_det[node]);
                problem.boundary_legal[a] = cost.exists && cost.w_int <= t_int;
                problem.boundary_w[a] = cost.w_int;
                // The §C adjacency carries both directions of every `H` edge and the local index of
                // a member is its position in this block, so one walk per member fills the pairing
                // table — the same legality test the decode path makes off `BallPrune`.
                for (uint32_t e = components.adj_offsets[node]; e < components.adj_offsets[node + 1]; e++) {
                    uint32_t b = components.local_index[components.adj_target[e]];
                    problem.pair_legal[a][b] = true;
                    problem.pair_w[a][b] = (pm::cumulative_time_int)components.adj_weight[e];
                }
            }
            verdict = resolve_small_component(problem).feasible ? SmallVerdict::COMMIT : SmallVerdict::RESIDUAL;
        }

        if (verdict == SmallVerdict::COMMIT) {
            stats.defects_committed_trivially += (int)size;
        } else if (verdict == SmallVerdict::RESIDUAL) {
            stats.defects_residual_trivially += (int)size;
        } else {
            stats.defects_to_solver += (int)size;
        }
    }

    // The three classifications partition `H`'s nodes; the aggregate's fractions are taken over
    // their sum, so a gap here would silently rescale every one of them.
    assert(
        stats.defects_committed_trivially + stats.defects_residual_trivially + stats.defects_to_solver ==
            (int)num_nodes &&
        "the resolver's classification did not partition H's defects");
}

void BallDecoder::compute_seeded_detection_events(const std::vector<uint64_t>& dets, std::vector<uint64_t>& out) const {
    const pm::MatchingGraph& graph = g_mwpm.flooder.graph;
    for (uint64_t det : dets) {
        if (det >= graph.nodes.size())
            throw std::invalid_argument(
                "The detection event with index " + std::to_string(det) +
                " does not correspond to a node in the graph, which only has " + std::to_string(graph.nodes.size()) +
                " nodes.");
    }

    // M1's negative-weight preamble, run on `G` *first*: a negative-weight edge implies a detection
    // event on both endpoints, and a shot's own event on such a node cancels it. `H` is then built
    // from the surviving set and never sees a negative weight (§M2.1).
    //
    // §B — on an all-positive DEM there is nothing to cancel against and the offsets are zero, so
    // the preamble is a no-op and the shot's own events go through directly. What is skipped is the
    // seeding alone; the boundary-node filter below is a separate requirement of `begin_timeline`'s
    // node set and runs either way.
    if (config.skip_negative_weight_preamble_when_positive && !dem_has_negative_weights) {
        assert(
            g_mwpm.flooder.negative_weight_detection_events.empty() &&
            "the preamble was skipped on a graph that has negative-weight detection events");
        out.assign(dets.begin(), dets.end());
        std::sort(out.begin(), out.end());
    } else {
        sort_scratch.assign(dets.begin(), dets.end());
        std::sort(sort_scratch.begin(), sort_scratch.end());
        const std::vector<uint64_t>& negatives = g_mwpm.flooder.negative_weight_detection_events;
        out.clear();
        std::set_symmetric_difference(
            sort_scratch.begin(), sort_scratch.end(), negatives.begin(), negatives.end(), std::back_inserter(out));
    }

    // A user-graph boundary node cannot carry a detection event; `begin_timeline` drops it, so `H`
    // must too or the two front ends would disagree on the node set.
    const auto& is_boundary = graph.is_user_graph_boundary_node;
    if (!is_boundary.empty()) {
        out.erase(
            std::remove_if(
                out.begin(),
                out.end(),
                [&](uint64_t det) {
                    return det < is_boundary.size() && is_boundary[det];
                }),
            out.end());
    }
}

template <typename HarvestOnH>
Phase1Outcome BallDecoder::decode_impl(
    const std::vector<uint64_t>& dets, BallProfile* prof, bool allow_prune, const HarvestOnH& harvest_on_h) {
    HiResTimer total_timer;
    HiResTimer step;
    if (prof != nullptr) {
        prof->clear();
        total_timer.start();
    }

    compute_seeded_detection_events(dets, seeded_scratch);

    // The structural counters need a profile *and* the flag: they are cheap enough to leave in the
    // profiled path only for the parts that are register increments, and the parts that are not
    // would move `intersect_ns` under their own measurement. See `BallConfig`.
    bool structural = prof != nullptr && config.collect_structural_counters;

    BallGraphTiming timing;
    BallGraphCounts counts;
    build_ball_graph(
        tables,
        seeded_scratch,
        horizon,
        arena,
        config.mode,
        prof != nullptr ? &timing : nullptr,
        structural ? &counts : nullptr);
    const BallGraph& h = arena.graph;

    // ---- §A. Resolve the components of size <= k off the ball tables and hand the solver what is
    // left.
    //
    // serial pre-pass; excluded from this measurement by intent. It runs before the solve, in
    // series, on the critical path, and it is a real cost — it is outside every `step.start()`
    // below because this branch's reported latency is scoped to the solver and the harvest on the
    // size-`> k` graph, to be measured and optimised separately, and not because it is free or
    // concurrent. `allow_prune` is false on the verification entry points, whose whole job is to
    // solve the same `H` the oracle does.
    resolved_pairs.clear();
    pm::MatchingResult resolved(0, 0);
    SmallCommitCounts resolved_counts;
    bool pruning = allow_prune && config.prune_component_max_size > 0;
    const BallGraph& solver_h = pruning ? resolve_small_components(h, resolved, resolved_counts) : h;
    solved_graph = &solver_h;

    // §A.5's verdict, known here rather than only at the combine below: a RESIDUAL component decides
    // the branch on its own, whatever the solver goes on to do with the size-`> k` remainder.
    //
    // It has to be known *before* the harvest. The production entry points abandon an escalating
    // shot in O(1) instead of extracting from it (§M3.4), and they decide that off the status they
    // are handed at the harvest call site — so leaving the residual to flip `outcome.status` after
    // that call had the shot extract a Phase 1 result that §A.5 then discarded, which is the wasted
    // work §M3.3 X9 counts.
    bool resolver_forces_escalation = pruning && resolved_counts.residual > 0;

    // §A.4. With every component resolved off the solver — the common case at `p = 1e-3` — there is no
    // `Mwpm(H)` to build, no timeline to run and nothing to extract. That is where the reduction
    // comes from, so it is a skip of the whole stage rather than a solve over an empty node set.
    //
    // The skip is conditioned on the prune, not merely on the node count, so that at `k = 0` the
    // path is byte-identical to the un-pruned one down to the harvest counters: a defect-free shot
    // still builds its empty `Mwpm(H)` and still runs its empty solve and extraction, exactly as it
    // does today.
    bool solver_runs = !pruning || solver_h.num_nodes() != 0;

    BallMwpmCounts mwpm_counts;
    if (solver_runs) {
        if (prof != nullptr)
            step.start();
        h_mwpm.rebuild(tables, solver_h, arena, structural ? &mwpm_counts : nullptr);
        if (prof != nullptr)
            prof->mwpm_build_ns = step.elapsed_ns();
    }

    // The solved graph's nodes are `0..n-1` by construction, so its detection events are every node.
    h_dets_scratch.clear();
    h_dets_scratch.reserve(solver_h.num_nodes());
    for (size_t i = 0; i < solver_h.num_nodes(); i++)
        h_dets_scratch.push_back(i);

    // §M2.9.6 measurement 3. A process-wide counter, so it is read as a delta around the solve.
    uint64_t formations_before = pm::blossom_formation_stats.formations;
    harvester.collect_diagnostics = config.collect_harvest_diagnostics;
    harvester.use_legacy_enumeration = config.use_legacy_harvest_enumeration;

    Phase1Outcome outcome;
    if (!solver_runs) {
        // §A.4's skip. Nothing ran, so nothing is reported: `blossom_on_h_ns`, `dual_scan_ns` and
        // `harvest_ns` stay at the zero `prof->clear()` left them at, which is the honest reading —
        // the stage did not happen. Under §M7 there is no dual scan to certify the shot with, and
        // none is needed *for the components the resolver committed*: it settled those exactly, at
        // `H`'s own cutoffs, so every committed match is inside the horizon the certificate tests —
        // a pair at `d <= 2T`, a boundary match at `bcost <= T`.
        //
        // This is the provisional reading only. A RESIDUAL component also leaves the solver empty
        // without being settled, and §A.5 below withdraws the certificate when there was one.
        if (prof != nullptr && config.stock_on_h)
            prof->certified = 1;
    } else if (config.stock_on_h) {
        // §M7. Stock blossom on `H`, no horizon, run to completion; the certificate is then read
        // once off the terminal dual. `TimelineStatus::TRUNCATED` keeps its meaning as *the branch*
        // — "Phase 1 has no usable answer, escalate" — which is all M3–M6 downstream consume; on
        // this path the timeline was never truncated, and `HarvestResult` is empty rather than
        // partial.
        if (prof != nullptr)
            step.start();
        CertificateOutcome certificate = run_stock_and_certify(h_mwpm.mwpm, h_dets_scratch, horizon, prof != nullptr);
        if (prof != nullptr) {
            // The solve proper, with the certificate's own scan netted out, so the two stages are
            // additive and §M7.8's read can charge the scan separately from the blossom work.
            prof->blossom_on_h_ns = step.elapsed_ns() - certificate.dual_scan_ns;
            prof->dual_scan_ns = certificate.dual_scan_ns;
            prof->certified = certificate.certified() ? 1 : 0;
            prof->h_no_perfect_matching = certificate.status == CertificateStatus::NO_PERFECT_MATCHING ? 1 : 0;
            prof->max_dual_at_completion = certificate.max_dual;
        }
        outcome.status = certificate.certified() ? TimelineStatus::COMPLETE : TimelineStatus::TRUNCATED;

        // Debug invariant 3: the shot escalates **iff** `H` had no perfect matching or it completed
        // with `max_u Y(u) > T_int`, read off the certificate's own recomputed dual rather than off
        // a proxy such as "the residual is empty" or "the solve threw".
        assert(
            (outcome.status == TimelineStatus::TRUNCATED) ==
                (certificate.status == CertificateStatus::NO_PERFECT_MATCHING || certificate.max_dual > horizon) &&
            "invariant 3: the escalation trigger and the certificate disagree");
        // Debug invariant 4, the machine guard against §M7.0's "catching the throw is enough"
        // fallacy: a run that *completed* on `H` with a dual over `T` must escalate. Trivially true
        // where it is written, and written anyway, because a refactor that reintroduced the fallacy
        // would be silent everywhere else.
        assert(
            !(certificate.status == CertificateStatus::DUAL_EXCEEDS_HORIZON &&
              outcome.status == TimelineStatus::COMPLETE) &&
            "invariant 4: a completing-but-over-T H-matching escaped as if it were certified");
    } else {
        if (prof != nullptr)
            step.start();
        outcome.status =
            config.collect_harvest_diagnostics
                ? process_timeline_until_horizon_measured(h_mwpm.mwpm, h_dets_scratch, horizon, depth_model)
                : process_timeline_until_horizon(h_mwpm.mwpm, h_dets_scratch, horizon);
        if (prof != nullptr)
            prof->blossom_on_h_ns = step.elapsed_ns();
    }

    // Invariant 5's second half: no exposed-root-blossom routine may be entered on a certified shot.
    // §M1.4 is unreachable here — its precondition is a *surviving* tree root — and this is that
    // statement asked of the code that would have run it rather than of the argument for why it
    // cannot.
    uint64_t base_descents_before = harvester.counters.base_descents;

    if (solver_runs) {
        if (prof != nullptr)
            step.start();
        // The status the *shot* has, not the one the solver's own components have: a residual
        // elsewhere means this result is about to be discarded, so the production path abandons
        // rather than extracts. `outcome.status` itself is left for §A.5, which is where the two
        // are reconciled and where the `stock_on_h` certificate has already written its verdict.
        outcome.harvest = harvest_on_h(
            h_mwpm.mwpm, h_dets_scratch, resolver_forces_escalation ? TimelineStatus::TRUNCATED : outcome.status);
        if (prof != nullptr)
            prof->harvest_ns = step.elapsed_ns();
    }
    HarvestResult& result = outcome.harvest;

    assert(
        (!config.stock_on_h || outcome.status != TimelineStatus::COMPLETE ||
         harvester.counters.base_descents == base_descents_before) &&
        "invariant 5: a certified shot entered the exposed-root-blossom base descent");
    (void)base_descents_before;

    // Back to `G`'s detector ids, through whichever graph the solver ran on. `h_to_det` is strictly
    // ascending in both cases — §A.4's sub-graph is a subsequence of `H` — so a residual sorted
    // there stays sorted in `G`, and nothing downstream of harvest has to know either graph existed.
    for (uint64_t& defect : result.residual) {
        assert(defect < solver_h.num_nodes());
        defect = solver_h.h_to_det[defect];
    }

    // ---- §A.5. Combine, and decide the branch.
    //
    // `escalate <=> resolver residual OR (the solver ran and truncated)`. That is the current
    // trigger restated, not a new one: the pipeline escalates a shot iff some component fails to
    // resolve within `T`, and here that surfaces either as a resolver RESIDUAL (a component of size
    // `<= k` with no feasible matching at all) or as the solver's own status on the components it
    // was handed. Nothing was deferred to the solver except by being size `> k`, so it decides
    // those exactly as it did before.
    if (pruning) {
        if (resolver_forces_escalation) {
            outcome.status = TimelineStatus::TRUNCATED;
            // §M7. The certificate is a statement about the *whole* of `H`, and the resolver owns
            // part of it, so its verdict has to reach the profile too. A RESIDUAL component is
            // exactly "`H` has no perfect matching": the enumeration is exhaustive over
            // pairings-with-boundary-fill at `H`'s own cutoffs, so no feasible configuration means
            // the solver — handed that component — would have reported NO_PERFECT_MATCHING itself.
            //
            // Both branches above can reach here with `certified == 1`: the empty-solver branch set
            // it provisionally, and the stock branch reads only the components it was handed, which
            // certify on their own while another component is residual. Withdrawing it here is what
            // keeps invariant 3's `escalated <=> !certified` true through the prune.
            //
            // `max_dual_at_completion` is deliberately left where it is. It is undefined rather
            // than merely unmeasured when `H` cannot complete (§M7), and the escalating shot
            // discards Phase 1 in full anyway.
            if (prof != nullptr && config.stock_on_h) {
                prof->certified = 0;
                prof->h_no_perfect_matching = 1;
            }
        }
        if (outcome.status == TimelineStatus::TRUNCATED) {
            // The escalating shot re-decodes the whole raw syndrome on `G` and discards **all** of
            // Phase 1, the off-solver commits included — there is nothing to XOR and nothing to add
            // (§A.5, §M3.3 X3). Dropping them here is what keeps a caller from combining them by
            // accident.
            resolved_pairs.clear();
        } else {
            result.committed.obs_mask ^= resolved.obs_mask;
            result.committed.weight += resolved.weight;
            // §M3.4/§M4.2. `dual_sum_at_truncation` is a sum over the regions of `H`, and the
            // resolver removed regions from the solve, so the harvest's reduction now runs over a
            // strict subset and under-reports the shot's dual by exactly what the resolver settled.
            //
            // What it settled is that dual: each committed match is frozen and tight, so the two
            // halves of a pair matched across `d` contribute `d / 2` each and a boundary match at
            // `bcost` contributes `bcost` — a component's regions sum to its committed weight,
            // which is `resolved.weight` over all of them. That is the same tightness §M4.2
            // measures as a ratio of exactly 1.0000 on a completed shot, and dropping the term
            // would have degraded the certificate silently: `weight_out >= dual_sum` survives a
            // dual that is too *small*, so only §M3.3 X8's bit-identity against the un-pruned full
            // harvest catches it.
            result.dual_sum_at_truncation += resolved.weight;
            // The off-solver commits are matches like any other, so they belong in the commit
            // tallies the profile reports; a pair settled by the resolver is frozen at the horizon
            // in the same sense a solver-frozen one is.
            result.committed_pairs_frozen += resolved_counts.pairs;
            result.committed_boundary += resolved_counts.boundary;
        }
    }

    if (prof != nullptr) {
        prof->intersect_ns = timing.intersect_ns;
        prof->h_build_ns = timing.finalize_ns;
        prof->n_defects = (int)dets.size();
        prof->h_nodes = (int)h.num_nodes();
        prof->h_edges = (int)h.edges.size();
        prof->h_boundary_edges = (int)h.boundary_edges.size();
        // The run's `k`, and what the resolver actually settled off the solver this shot. Both are
        // labels on the latency numbers beside them, not latency numbers themselves: the resolve is
        // a serial pre-pass excluded from `blossom_on_h_ns` and `harvest_ns` by scope. `k` is
        // reported as the decoder was configured, so a `k > 0` run whose shot resolved nothing is
        // still distinguishable from a `k = 0` one.
        prof->prune_component_max_size = pruning ? config.prune_component_max_size : 0;
        prof->defects_resolved_small = resolved_counts.defects_resolved;
        // §M2 structural counters, left at zero unless they were collected, so that a profile
        // never reports a counter it did not measure. `hbld_edges_written` is the count taken at
        // the `push_back`s, not `h_edges + h_boundary_edges` read back off the vectors —
        // `build_ball_graph` asserts the two agree, which is where the "each undirected pair once"
        // claim is checked.
        if (structural) {
            prof->isect_scan_bytes = counts.isect_scan_bytes;
            prof->isect_hit_bytes = counts.isect_hit_bytes;
            prof->isect_scan_bytes_other_mode = counts.isect_scan_bytes_other_mode;
            prof->isect_mode_bitset = config.mode == BallGraphBuildMode::BITSET;
            prof->hbld_edges_written = (int)(counts.edges_written + counts.boundary_edges_written);
            prof->mwpm_init_node_elements = (int)mwpm_counts.node_records;
            prof->mwpm_init_edge_elements = (int)mwpm_counts.edge_records;
        }
        prof->shells_materialized = 1;
        prof->restarts = 0;
        prof->blossom_formations = (int)(pm::blossom_formation_stats.formations - formations_before);
        prof->harvest_enumerate_ns = result.enumerate_ns;
        prof->harvest_reduce_ns = result.reduce_ns;
        prof->harvest_base_descent_ns = result.base_descent_ns;
        prof->harvest_shatter_ns = result.shatter_ns;
        prof->max_blossom_nesting_depth = result.max_blossom_nesting_depth;
        prof->max_blossom_members = result.max_blossom_members;
        prof->max_exposed_blossom_depth = result.max_exposed_blossom_depth;
        prof->max_exposed_blossom_members = result.max_exposed_blossom_members;
        prof->matched_blossom_shatters = result.matched_blossom_shatters;
        prof->largest_tree_size = result.largest_tree_size;
        prof->harvest_dependent_depth = result.harvest_dependent_depth;
        // §M2.9.6 measurement 4 is a property of the *truncated* timeline loop, which is the only
        // one instrumented; §M7's front end runs stock's loop untouched, which is the point, so it
        // reports no solve depth rather than a stale one from the previous shot.
        bool measured_solve = config.collect_harvest_diagnostics && !config.stock_on_h;
        prof->solve_dependent_depth = measured_solve ? depth_model.depth : 0;
        prof->solve_events = measured_solve ? depth_model.events : 0;
        if (h.num_nodes() != 0)
            prof->mean_degree = 2.0 * (double)h.edges.size() / (double)h.num_nodes();
        // `max_degree` is read off the matching graph the solve was built on, so under §A it is the
        // largest degree the *solver* saw rather than `H`'s. That is the quantity the number is
        // used for; `H`'s own degree distribution is §C's `degree_hist`, which is taken over the
        // whole graph whether or not anything was pruned.
        if (solver_runs) {
            for (size_t i = 0; i < solver_h.num_nodes(); i++)
                prof->max_degree = std::max(prof->max_degree, (int)h_mwpm.mwpm.flooder.graph.nodes[i].neighbors.size());
        }
        prof->total_ns = total_timer.elapsed_ns();
    }
    return outcome;
}

void BallDecoder::map_match_edges_to_committed_pairs(std::vector<CommittedPair>& committed_pairs) const {
    // The graph the solve ran on: the whole of `H`, or §A.4's sub-`H`. The solver's node indices are
    // indices into *that* graph, and so is the edge list the ball entry is looked up in.
    const BallGraph& h = solved_graph != nullptr ? *solved_graph : arena.graph;
    const pm::DetectorNode* base = h_mwpm.mwpm.flooder.graph.nodes.data();
    committed_pairs.clear();
    committed_pairs.reserve(match_edge_scratch.size() + resolved_pairs.size());
    for (const pm::CompressedEdge& edge : match_edge_scratch) {
        size_t i = (size_t)(edge.loc_from - base);
        int64_t from = (int64_t)h.h_to_det[i];
        if (edge.loc_to == nullptr) {
            committed_pairs.push_back(CommittedPair{from, -1, CommittedPair::NO_BALL_ENTRY});
            continue;
        }
        size_t j = (size_t)(edge.loc_to - base);
        int64_t to = (int64_t)h.h_to_det[j];
        // Every committed pair is an edge of `H` — a region only ever meets another region across
        // one — so the ball entry behind it is in `h.edges`, which is sorted by `(i, j)` with
        // `i < j`. One binary search; no map, no per-shot allocation.
        uint32_t lo = (uint32_t)std::min(i, j);
        uint32_t hi = (uint32_t)std::max(i, j);
        auto it = std::lower_bound(
            h.edges.begin(),
            h.edges.end(),
            std::pair<uint32_t, uint32_t>{lo, hi},
            [](const BallGraphEdge& e, const std::pair<uint32_t, uint32_t>& key) {
                return e.i != key.first ? e.i < key.first : e.j < key.second;
            });
        assert(it != h.edges.end() && it->i == lo && it->j == hi && "a committed pair that is not an edge of H");
        committed_pairs.push_back(CommittedPair{from, to, it->entry});
    }
    // §A's off-solver commits, already in `G`'s detector ids and already carrying the ball entry
    // the mask and weight came from. Empty unless the prune ran and the shot completed.
    committed_pairs.insert(committed_pairs.end(), resolved_pairs.begin(), resolved_pairs.end());
    sort_pairs(committed_pairs);
}

namespace {

/// Why the two harvest-flavoured entry points are closed on the stock path. See the header.
[[noreturn]] void reject_harvest_entry_point() {
    throw std::invalid_argument(
        "decode_phase1/decode_phase1_to_match_edges return a HarvestResult, and their callers read the escalation "
        "decision off an empty residual. Under the §M7 certificate that reading is wrong exactly where it matters: "
        "a shot that completes on H with a suboptimal matching and max_u Y(u) > T has an empty residual and must "
        "still escalate. Use decode_phase1_production*(), which carries the certificate's own status.");
}

}  // namespace

HarvestResult BallDecoder::decode_phase1(const std::vector<uint64_t>& dets, BallProfile* prof) {
    if (config.stock_on_h)
        reject_harvest_entry_point();
    // No prune (§A): this is the verification path, and the oracle it is compared against solves
    // the whole of `H`.
    Phase1Outcome outcome =
        decode_impl(dets, prof, false, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus) {
            return harvester.harvest_to_obs(mwpm, h_dets);
        });
    if (config.verify_against_g)
        verify_level1(dets, outcome.harvest, nullptr, prof);
    return outcome.harvest;
}

HarvestResult BallDecoder::decode_phase1_to_match_edges(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof) {
    if (config.stock_on_h)
        reject_harvest_entry_point();
    match_edge_scratch.clear();
    Phase1Outcome outcome =
        decode_impl(dets, prof, false, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus) {
            return harvester.harvest_to_match_edges(mwpm, h_dets, match_edge_scratch);
        });
    map_match_edges_to_committed_pairs(committed_pairs);

    if (config.verify_against_g)
        verify_level1(dets, outcome.harvest, &committed_pairs, prof);
    return outcome.harvest;
}

Phase1Outcome BallDecoder::decode_phase1_production(const std::vector<uint64_t>& dets, BallProfile* prof) {
    return decode_impl(
        dets, prof, true, [this, prof](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus status) {
            if (status == TimelineStatus::TRUNCATED) {
                HiResTimer abandon;
                if (prof != nullptr)
                    abandon.start();
                abandon_shot(mwpm);
                if (prof != nullptr)
                    prof->abandon_ns = abandon.elapsed_ns();
                return HarvestResult();
            }
            return harvester.extract_only_to_obs(mwpm, h_dets);
        });
}

Phase1Outcome BallDecoder::decode_phase1_production_to_match_edges(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof) {
    match_edge_scratch.clear();
    Phase1Outcome outcome = decode_impl(
        dets, prof, true, [this, prof](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus status) {
            if (status == TimelineStatus::TRUNCATED) {
                HiResTimer abandon;
                if (prof != nullptr)
                    abandon.start();
                abandon_shot(mwpm);
                if (prof != nullptr)
                    prof->abandon_ns = abandon.elapsed_ns();
                return HarvestResult();
            }
            return harvester.extract_only_to_match_edges(mwpm, h_dets, match_edge_scratch);
        });
    // On an escalating shot `match_edge_scratch` is empty and this yields no pairs, which is the
    // right answer: the caller is about to discard Phase 1 entirely.
    map_match_edges_to_committed_pairs(committed_pairs);
    return outcome;
}

bool BallDecoder::truncated_scheme_escalates(const std::vector<uint64_t>& dets) {
    // The same `H`, built from the same tables at the same `T` filter — the two schemes differ only
    // in the front end, so replaying the decision means replaying the solve, not rebuilding the
    // problem differently.
    compute_seeded_detection_events(dets, seeded_scratch);
    build_ball_graph(tables, seeded_scratch, horizon, arena, config.mode, nullptr, nullptr);
    const BallGraph& h = arena.graph;
    // The landed scheme is what is being replayed, so no §A prune: it solves the whole of `H`. The
    // pointer the last real decode left behind is dropped with it — the graph it named has just been
    // rebuilt underneath it, and nothing here produces pairs to map through it.
    solved_graph = nullptr;
    h_mwpm.rebuild(tables, h, arena, nullptr);

    h_dets_scratch.clear();
    h_dets_scratch.reserve(h.num_nodes());
    for (size_t i = 0; i < h.num_nodes(); i++)
        h_dets_scratch.push_back(i);

    TimelineStatus status = process_timeline_until_horizon(h_mwpm.mwpm, h_dets_scratch, horizon);
    // Nothing is harvested and nothing is read: only the *decision* is wanted. `abandon_shot` is
    // the teardown that works from either outcome, and it is what an escalating shot pays anyway.
    abandon_shot(h_mwpm.mwpm);
    return status == TimelineStatus::TRUNCATED;
}

HarvestResult BallDecoder::reference_phase1_on_g(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>* committed_pairs) {
    if (committed_pairs == nullptr) {
        process_timeline_until_horizon(g_mwpm, dets, horizon);
        return harvester.harvest_to_obs(g_mwpm, dets);
    }
    std::vector<pm::CompressedEdge> match_edges;
    process_timeline_until_horizon(g_mwpm, dets, horizon);
    HarvestResult result = harvester.harvest_to_match_edges(g_mwpm, dets, match_edges);
    const pm::DetectorNode* base = g_mwpm.flooder.graph.nodes.data();
    committed_pairs->clear();
    committed_pairs->reserve(match_edges.size());
    for (const pm::CompressedEdge& edge : match_edges) {
        int64_t from = (int64_t)(edge.loc_from - base);
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)(edge.loc_to - base);
        committed_pairs->push_back(CommittedPair{from, to});
    }
    sort_pairs(*committed_pairs);
    return result;
}

void BallDecoder::verify_level1(
    const std::vector<uint64_t>& dets,
    const HarvestResult& actual,
    const std::vector<CommittedPair>* actual_pairs,
    BallProfile* prof) {
    HiResTimer timer;
    timer.start();
    std::vector<CommittedPair> expected_pairs;
    HarvestResult expected = reference_phase1_on_g(dets, actual_pairs != nullptr ? &expected_pairs : nullptr);
    if (prof != nullptr)
        prof->g_reference_ns = timer.elapsed_ns();

    // §M2.6 level 1. Everything the dual solution determines is compared with no tolerance, and a
    // difference throws where it happens — this is debug invariant 10 and a release blocker.
    //
    // The design's list also asks for the residual *set* and the committed *pair set* to match
    // outright, reasoning that "degenerate path choice cannot move them". Path choice cannot, but
    // the choice *among optimal primal solutions* can, for a reason §M2.0 does not cover: in `G` a
    // growing region's flood is blocked by its neighbours' territory, so `G` never observes some
    // tight collisions that `H` — where every pair within `2T` is a direct edge — does observe. Both
    // are correct blossom implementations of the same metric; they agree on the dual solution and on
    // the optimum, not necessarily on which optimum they land on when tight events coincide.
    //
    // Measured over 4800 shot-decodes across six corpora and four horizons: the dual sum, the
    // committed weight, `num_trees` and the observable bytes matched on every one; only the pairing
    // (0–25% of shots, rising with defect density) and, at `T <= 1.0` on the densest corpora only,
    // the residual choice moved. Both residuals always satisfied `Y(u) == T`.
    //
    // `committed_boundary` belongs with those: it is a property of the chosen pairing, not of the
    // dual solution. When a tight collision lets `H` pair two defects that `G` matched to the
    // boundary separately — same total weight, same `num_trees`, same optimum — the boundary count
    // moves with the pairing. It matched on all 4800 of those shot-decodes, but nothing forces it
    // to, so it is counted rather than thrown on.
    //
    // Those three are therefore counted, not thrown on; everything else is still fatal.
    auto fail = [&](const std::string& what, const std::string& expected_text, const std::string& actual_text) {
        throw std::logic_error(
            "Ball-graph decode diverged from M1 on G at §M2.6 level 1 (" + what + "): G gave " + expected_text +
            ", H gave " + actual_text + ". This is debug invariant 10 and a release blocker.");
    };
    if (actual.dual_sum_at_truncation != expected.dual_sum_at_truncation)
        fail(
            "dual_sum_at_truncation",
            std::to_string(expected.dual_sum_at_truncation),
            std::to_string(actual.dual_sum_at_truncation));
    if (actual.committed.weight != expected.committed.weight)
        fail("committed weight", std::to_string(expected.committed.weight), std::to_string(actual.committed.weight));
    if (actual.num_trees != expected.num_trees)
        fail("num_trees", std::to_string(expected.num_trees), std::to_string(actual.num_trees));
    if (actual.residual.size() != expected.residual.size())
        fail("residual size", describe(expected.residual), describe(actual.residual));
    // The separation invariant, which is what says a differently-chosen residual is still a valid
    // exposed set. If this fires, the Phase-2 error bound is void.
    for (pm::total_weight_int y : actual.residual_dual_sum) {
        if (y != horizon)
            fail("separation invariant Y(u) == T", std::to_string(horizon), std::to_string(y));
    }

    bool residual_tie = actual.residual != expected.residual;
    bool boundary_tie = actual.committed_boundary != expected.committed_boundary;
    bool pairing_tie = false;
    if (actual_pairs != nullptr) {
        // The pair *count* is not independent of the boundary count: with the same number of
        // matched defects, `pairs == (matched + committed_boundary) / 2`, so the two differ on
        // exactly the same shots. Comparing it separately would reintroduce the check above under
        // another name, so a difference here is the same tie.
        //
        // Debug invariant 3: committed and residual partition the shot's detection events, each
        // classified exactly once. The committed support alone is *not* comparable across the two
        // front ends — a differently resolved residual tie moves it by the complementary defect —
        // but the partition is invariant, and it is what says nothing was dropped or double-counted.
        std::vector<uint64_t> partition(actual.residual);
        for (const CommittedPair& pair : *actual_pairs) {
            partition.push_back((uint64_t)pair.from);
            if (pair.to >= 0)
                partition.push_back((uint64_t)pair.to);
        }
        std::sort(partition.begin(), partition.end());
        if (partition != seeded_scratch)
            fail("syndrome partition", describe(seeded_scratch), describe(partition));
        pairing_tie = !std::equal(
            actual_pairs->begin(),
            actual_pairs->end(),
            expected_pairs.begin(),
            expected_pairs.end(),
            [](const CommittedPair& a, const CommittedPair& b) {
                return a.from == b.from && a.to == b.to;
            });
    }
    if (prof != nullptr) {
        prof->residual_ties += residual_tie ? 1 : 0;
        prof->pairing_ties += pairing_tie ? 1 : 0;
        prof->boundary_ties += boundary_tie ? 1 : 0;
    }
}

}  // namespace spec_matching
}  // namespace pm
