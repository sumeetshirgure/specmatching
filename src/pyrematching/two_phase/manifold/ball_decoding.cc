// Copyright 2026 PyReMatching contributors
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

#include "pyrematching/two_phase/manifold/ball_decoding.h"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <sstream>

#include "pyrematching/two_phase/manifold/ball_serialize.h"

namespace pm {
namespace two_phase {

namespace {

std::string describe(const std::vector<uint64_t>& values) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < values.size(); i++)
        out << (i ? ", " : "") << values[i];
    out << "]";
    return out.str();
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

    h_mwpm.configure(graph.num_observables, graph.normalising_constant);
    arena.reset_for_graph(graph.nodes.size());
}

void BallDecoder::save_ball_artifact(const std::string& path) const {
    save_ball_tables(tables, path);
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
    sort_scratch.assign(dets.begin(), dets.end());
    std::sort(sort_scratch.begin(), sort_scratch.end());
    const std::vector<uint64_t>& negatives = g_mwpm.flooder.negative_weight_detection_events;
    out.clear();
    std::set_symmetric_difference(
        sort_scratch.begin(), sort_scratch.end(), negatives.begin(), negatives.end(), std::back_inserter(out));

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
HarvestResult BallDecoder::decode_impl(
    const std::vector<uint64_t>& dets, BallProfile* prof, const HarvestOnH& harvest_on_h) {
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

    BallMwpmCounts mwpm_counts;
    if (prof != nullptr)
        step.start();
    h_mwpm.rebuild(tables, h, arena, structural ? &mwpm_counts : nullptr);
    if (prof != nullptr)
        prof->mwpm_build_ns = step.elapsed_ns();

    // `H`'s nodes are `0..n-1` by construction, so its detection events are every node.
    h_dets_scratch.clear();
    h_dets_scratch.reserve(h.num_nodes());
    for (size_t i = 0; i < h.num_nodes(); i++)
        h_dets_scratch.push_back(i);

    // §M2.9.6 measurement 3. A process-wide counter, so it is read as a delta around the solve.
    uint64_t formations_before = pm::blossom_formation_stats.formations;
    harvester.collect_diagnostics = config.collect_harvest_diagnostics;
    harvester.use_legacy_enumeration = config.use_legacy_harvest_enumeration;

    if (prof != nullptr)
        step.start();
    TimelineStatus status = config.collect_harvest_diagnostics
                                ? process_timeline_until_horizon_measured(
                                      h_mwpm.mwpm, h_dets_scratch, horizon, depth_model)
                                : process_timeline_until_horizon(h_mwpm.mwpm, h_dets_scratch, horizon);
    (void)status;
    if (prof != nullptr)
        prof->blossom_on_h_ns = step.elapsed_ns();

    if (prof != nullptr)
        step.start();
    HarvestResult result = harvest_on_h(h_mwpm.mwpm, h_dets_scratch);
    if (prof != nullptr)
        prof->harvest_ns = step.elapsed_ns();

    // Back to `G`'s detector ids. `h_to_det` is strictly ascending, so a residual sorted in `H`
    // stays sorted in `G` — nothing downstream of harvest has to know `H` existed.
    for (uint64_t& defect : result.residual) {
        assert(defect < h.num_nodes());
        defect = h.h_to_det[defect];
    }

    if (prof != nullptr) {
        prof->intersect_ns = timing.intersect_ns;
        prof->h_build_ns = timing.finalize_ns;
        prof->n_defects = (int)dets.size();
        prof->h_nodes = (int)h.num_nodes();
        prof->h_edges = (int)h.edges.size();
        prof->h_boundary_edges = (int)h.boundary_edges.size();
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
        prof->solve_dependent_depth = config.collect_harvest_diagnostics ? depth_model.depth : 0;
        prof->solve_events = config.collect_harvest_diagnostics ? depth_model.events : 0;
        if (h.num_nodes() != 0) {
            prof->mean_degree = 2.0 * (double)h.edges.size() / (double)h.num_nodes();
            for (size_t i = 0; i < h.num_nodes(); i++) {
                prof->max_degree = std::max(prof->max_degree, (int)h_mwpm.mwpm.flooder.graph.nodes[i].neighbors.size());
            }
        }
        prof->total_ns = total_timer.elapsed_ns();
    }
    return result;
}

HarvestResult BallDecoder::decode_phase1(const std::vector<uint64_t>& dets, BallProfile* prof) {
    HarvestResult result = decode_impl(dets, prof, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets) {
        return harvester.harvest_to_obs(mwpm, h_dets);
    });
    if (config.verify_against_g)
        verify_level1(dets, result, nullptr, prof);
    return result;
}

HarvestResult BallDecoder::decode_phase1_to_match_edges(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof) {
    match_edge_scratch.clear();
    HarvestResult result = decode_impl(dets, prof, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets) {
        return harvester.harvest_to_match_edges(mwpm, h_dets, match_edge_scratch);
    });

    const BallGraph& h = arena.graph;
    const pm::DetectorNode* base = h_mwpm.mwpm.flooder.graph.nodes.data();
    committed_pairs.clear();
    committed_pairs.reserve(match_edge_scratch.size());
    for (const pm::CompressedEdge& edge : match_edge_scratch) {
        int64_t from = (int64_t)h.h_to_det[(size_t)(edge.loc_from - base)];
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)h.h_to_det[(size_t)(edge.loc_to - base)];
        committed_pairs.push_back(CommittedPair{from, to});
    }
    sort_pairs(committed_pairs);

    if (config.verify_against_g)
        verify_level1(dets, result, &committed_pairs, prof);
    return result;
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

}  // namespace two_phase
}  // namespace pm
