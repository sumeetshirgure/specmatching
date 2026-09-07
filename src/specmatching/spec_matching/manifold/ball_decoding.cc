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

/// §2's per-component harvests, folded into the shot's one `HarvestResult`.
///
/// Every field here is either a XOR, a sum or a max over the components, which is what §1's
/// independence property says it has to be: the components' event sequences are disjoint, so the
/// shot's committed observable is the XOR of theirs, its weight and dual the sums, and its residual
/// the union. `part` is moved from, so the caller must not read it afterwards.
void accumulate_component_harvest(HarvestResult& into, HarvestResult& part) {
    into.committed.obs_mask ^= part.committed.obs_mask;
    into.committed.weight += part.committed.weight;

    into.residual.insert(into.residual.end(), part.residual.begin(), part.residual.end());
    into.residual_dual_sum.insert(
        into.residual_dual_sum.end(), part.residual_dual_sum.begin(), part.residual_dual_sum.end());

    into.committed_pairs_frozen += part.committed_pairs_frozen;
    into.committed_pairs_tree += part.committed_pairs_tree;
    into.committed_pairs_blossom_cycle += part.committed_pairs_blossom_cycle;
    into.committed_boundary += part.committed_boundary;
    into.num_trees += part.num_trees;
    into.exposed_root_blossoms += part.exposed_root_blossoms;

    into.dual_sum_at_truncation += part.dual_sum_at_truncation;
    into.max_region_dual = std::max(into.max_region_dual, part.max_region_dual);

    // Maxima, because they describe the worst structure anywhere in the shot — the same reading
    // they had when one solve saw all of `H`.
    into.largest_tree_size = std::max(into.largest_tree_size, part.largest_tree_size);
    into.max_blossom_nesting_depth = std::max(into.max_blossom_nesting_depth, part.max_blossom_nesting_depth);
    into.max_blossom_members = std::max(into.max_blossom_members, part.max_blossom_members);
    into.max_exposed_blossom_depth = std::max(into.max_exposed_blossom_depth, part.max_exposed_blossom_depth);
    into.max_exposed_blossom_members = std::max(into.max_exposed_blossom_members, part.max_exposed_blossom_members);
    into.matched_blossom_shatters += part.matched_blossom_shatters;
    into.harvest_dependent_depth = std::max(into.harvest_dependent_depth, part.harvest_dependent_depth);

    into.enumerate_ns += part.enumerate_ns;
    into.reduce_ns += part.reduce_ns;
    into.base_descent_ns += part.base_descent_ns;
    into.shatter_ns += part.shatter_ns;
}

/// `HarvestResult::residual` is contracted to be sorted ascending, and the union of the components'
/// residuals is not: each component's is sorted within itself, but the blocks interleave. One sort
/// of the `(defect, Y(u))` pairs restores it, which is also what makes the §2.5 set comparison
/// against the monolithic residual a plain `==`.
void sort_residual(HarvestResult& result, std::vector<std::pair<uint64_t, pm::total_weight_int>>& scratch) {
    if (result.residual.size() < 2)
        return;
    scratch.clear();
    scratch.reserve(result.residual.size());
    for (size_t i = 0; i < result.residual.size(); i++)
        scratch.emplace_back(result.residual[i], result.residual_dual_sum[i]);
    std::sort(scratch.begin(), scratch.end());
    for (size_t i = 0; i < scratch.size(); i++) {
        result.residual[i] = scratch[i].first;
        result.residual_dual_sum[i] = scratch[i].second;
    }
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
    // §2.5's monolithic instance is configured either way — it costs one empty `pm::Mwpm` — but it
    // only ever grows a node pool if the flag actually turns the check on.
    verify_mwpm.configure(graph.num_observables, graph.normalising_constant);
    arena.reset_for_graph(graph.nodes.size());
}

void BallDecoder::set_horizon(double T) {
    if (T > config.ball.T_max)
        throw std::invalid_argument(
            "BallConfig::T exceeds BallParams::T_max; the compiled tables do not cover that horizon.");
    horizon_int next = to_time_units(T, g_mwpm.flooder.graph.normalising_constant);
    if (next > tables.t_max_int)
        throw std::invalid_argument("BallConfig::T exceeds the compiled T_max in time units.");
    config.T = T;
    horizon = next;
}

void BallDecoder::save_ball_artifact(const std::string& path) const {
    save_ball_tables(tables, path);
}

void BallDecoder::analyze_last_shot_components(
    BallProfile& prof, ComponentHistograms& histograms, ComponentStatusTable& size_x_status) {
    // No timer is started anywhere in this function, and none may be added: §3.5.2 is untimed by
    // construction and the caller has already closed the shot's window.
    ComponentStats& stats = prof.components;
    stats = ComponentStats();
    histograms.clear();
    size_x_status.clear();

    const BallGraph& h = arena.graph;
    BallComponents& components = arena.components;
    analyze_ball_components(h, components);

    // The decode's own statuses, for the joint table. Both decompositions are union-find over the
    // same edge set enumerated in ascending root order, so component `c` is the same component in
    // both — asserted rather than assumed, because a joint table joined on the wrong key is a table
    // that looks right and says nothing.
    const BallComponentSplit& split = arena.split;
    assert(split.roots.size() == components.roots.size() && "the two decompositions of H disagree on the components");
    assert(
        std::equal(split.roots.begin(), split.roots.end(), components.roots.begin()) &&
        "the two decompositions of H disagree on the component roots");

    // The horizon in the decoder's own stored integer units. Read, never re-derived: converting
    // `config.T` a second time here would be the §0 unit-rule trap.
    const pm::cumulative_time_int t_int = horizon;
    uint32_t num_nodes = (uint32_t)h.num_nodes();

    stats.measured = 1;
    stats.num_components = (int)components.num_components();
    stats.component_defects = (int)num_nodes;
    stats.components_truncated = (int)split.num_truncated();

    // ---- Per-edge and per-defect distributions. Every `H` edge weight is an exact `d_G` between two
    // defects, which is the path-length distribution §C.2 asks for. These two are the latency
    // profiler's accumulators — they are keyed by edge and by defect, not by component, which is why
    // they survived the cut to size-only component statistics. `sparse_graph_stats` does not write
    // them (§3.5.2, "Not collected").
    for (const BallGraphEdge& edge : h.edges) {
        histograms.edge_weight_hist[ComponentHistograms::weight_bin(
            (pm::cumulative_time_int)edge.w_int, t_int, histograms.edge_weight_hist.size())]++;
    }
    for (uint32_t i = 0; i < num_nodes; i++) {
        BoundaryCost cost = boundary_cost(tables, h.h_to_det[i]);
        if (cost.exists) {
            histograms.bcost_hist[ComponentHistograms::weight_bin(
                cost.w_int, t_int, histograms.bcost_hist.size())]++;
        }
    }

    // ---- Per component, in ascending root order (§0 determinism). Size and status, and nothing
    // else: no adjacency is built, no member block is walked, and no diameter is computed.
    for (size_t c = 0; c < components.num_components(); c++) {
        uint32_t size = components.sizes[c];
        size_t status = split.status[c] == BallComponentSplit::COMPONENT_TRUNCATED ? ComponentStatusTable::TRUNCATED
                                                                                  : ComponentStatusTable::COMPLETE;

        histograms.add_component_size(size);
        size_x_status.add(size, status);
        stats.largest_component_size = std::max(stats.largest_component_size, (int)size);

        // The size classes, which are structural: "trivial" is size <= 2 whatever the solve does
        // with it, so the series stays comparable across runs and across branches.
        if (size == 1) {
            stats.num_singleton_components++;
        } else if (size == 2) {
            stats.num_pair_components++;
        } else {
            stats.num_components_size_ge3++;
        }
        if (size <= 2) {
            stats.num_trivial_components++;
            stats.defects_in_trivial_components += (int)size;
        }
    }

    // The components partition `H`'s nodes; every per-shot fraction below is taken over that sum, so
    // a gap here would silently rescale all of them.
    assert(
        stats.num_singleton_components + 2 * stats.num_pair_components <= (int)num_nodes &&
        "the component size classes overran H's node count");
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

void BallDecoder::drain_component_match_edges() {
    // The sub-`H` the component instance was just built on. The instance's node indices are indices
    // into *that* graph, and so is the edge list the ball entry is looked up in — which is why this
    // has to run before the next component rebuilds the instance underneath these pointers.
    const BallGraph& sub = arena.split.sub;
    const pm::DetectorNode* base = h_mwpm.mwpm.flooder.graph.nodes.data();
    for (const pm::CompressedEdge& edge : match_edge_scratch) {
        size_t i = (size_t)(edge.loc_from - base);
        int64_t from = (int64_t)sub.h_to_det[i];
        if (edge.loc_to == nullptr) {
            committed_pair_scratch.push_back(CommittedPair{from, -1, CommittedPair::NO_BALL_ENTRY});
            continue;
        }
        size_t j = (size_t)(edge.loc_to - base);
        int64_t to = (int64_t)sub.h_to_det[j];
        // Every committed pair is an edge of the component's sub-`H` — a region only ever meets
        // another region across one, and §1 says it cannot meet one outside its component at all —
        // so the ball entry behind it is in `sub.edges`, which is sorted by `(i, j)` with `i < j`.
        // One binary search; no map, no per-shot allocation.
        uint32_t lo = (uint32_t)std::min(i, j);
        uint32_t hi = (uint32_t)std::max(i, j);
        auto it = std::lower_bound(
            sub.edges.begin(),
            sub.edges.end(),
            std::pair<uint32_t, uint32_t>{lo, hi},
            [](const BallGraphEdge& e, const std::pair<uint32_t, uint32_t>& key) {
                return e.i != key.first ? e.i < key.first : e.j < key.second;
            });
        assert(it != sub.edges.end() && it->i == lo && it->j == hi && "a committed pair that is not an edge of H[C]");
        committed_pair_scratch.push_back(CommittedPair{from, to, it->entry});
    }
    match_edge_scratch.clear();
}

template <typename HarvestComponent>
TimelineStatus BallDecoder::solve_component(
    BallProfile* prof, const HarvestComponent& harvest_component, HarvestResult& harvest) {
    const BallGraph& sub = arena.split.sub;
    HiResTimer step;

    if (prof != nullptr)
        step.start();
    BallMwpmCounts mwpm_counts;
    bool structural = prof != nullptr && config.collect_structural_counters;
    h_mwpm.rebuild(tables, sub, arena, structural ? &mwpm_counts : nullptr);
    if (prof != nullptr) {
        prof->mwpm_build_ns += step.elapsed_ns();
        if (structural) {
            prof->mwpm_init_node_elements += (int)mwpm_counts.node_records;
            prof->mwpm_init_edge_elements += (int)mwpm_counts.edge_records;
        }
    }

    // The component's nodes are `0..s-1` by construction, so its detection events are every node.
    h_dets_scratch.clear();
    h_dets_scratch.reserve(sub.num_nodes());
    for (size_t i = 0; i < sub.num_nodes(); i++)
        h_dets_scratch.push_back(i);

    TimelineStatus status = TimelineStatus::COMPLETE;
    if (config.stock_on_h) {
        // §M7's front end, asked of the component rather than of `H`. The certificate composes over
        // the split for the same reason the timeline does: `max_u Y(u)` over `H` is the max over the
        // components of their own terminal duals, and `H` has a perfect matching iff every component
        // does. Both follow from §1 — the regions of two components never interact before `T`.
        if (prof != nullptr)
            step.start();
        CertificateOutcome certificate = run_stock_and_certify(h_mwpm.mwpm, h_dets_scratch, horizon, prof != nullptr);
        if (prof != nullptr) {
            // The solve proper, with the certificate's own scan netted out, so the two stages are
            // additive and §M7.8's read can charge the scan separately from the blossom work.
            prof->blossom_on_h_ns += step.elapsed_ns() - certificate.dual_scan_ns;
            prof->dual_scan_ns += certificate.dual_scan_ns;
            prof->h_no_perfect_matching |= certificate.status == CertificateStatus::NO_PERFECT_MATCHING ? 1 : 0;
            prof->max_dual_at_completion = std::max(prof->max_dual_at_completion, certificate.max_dual);
        }
        status = certificate.certified() ? TimelineStatus::COMPLETE : TimelineStatus::TRUNCATED;

        // Debug invariant 3, per component: it escalates **iff** this component had no perfect
        // matching or completed with `max_u Y(u) > T_int`, read off the certificate's own
        // recomputed dual rather than off a proxy such as "the residual is empty".
        assert(
            (status == TimelineStatus::TRUNCATED) ==
                (certificate.status == CertificateStatus::NO_PERFECT_MATCHING || certificate.max_dual > horizon) &&
            "invariant 3: the escalation trigger and the certificate disagree");
        // Debug invariant 4, the machine guard against §M7.0's "catching the throw is enough"
        // fallacy: a component that *completed* with a dual over `T` must escalate.
        assert(
            !(certificate.status == CertificateStatus::DUAL_EXCEEDS_HORIZON && status == TimelineStatus::COMPLETE) &&
            "invariant 4: a completing-but-over-T H-matching escaped as if it were certified");
    } else {
        if (prof != nullptr)
            step.start();
        status = config.collect_harvest_diagnostics
                     ? process_timeline_until_horizon_measured(h_mwpm.mwpm, h_dets_scratch, horizon, depth_model)
                     : process_timeline_until_horizon(h_mwpm.mwpm, h_dets_scratch, horizon);
        if (prof != nullptr) {
            prof->blossom_on_h_ns += step.elapsed_ns();
            if (config.collect_harvest_diagnostics) {
                // §M2.9.6 measurement 4 is a serial depth, so it is the deepest chain anywhere in
                // the shot; the event count is the total.
                prof->solve_dependent_depth = std::max(prof->solve_dependent_depth, depth_model.depth);
                prof->solve_events += depth_model.events;
            }
        }
    }

    // Invariant 5's second half: no exposed-root-blossom routine may be entered on a certified
    // component. §M1.4 is unreachable there — its precondition is a *surviving* tree root.
    uint64_t base_descents_before = harvester.counters.base_descents;

    if (prof != nullptr)
        step.start();
    match_edge_scratch.clear();
    harvest = harvest_component(h_mwpm.mwpm, h_dets_scratch, status);
    if (prof != nullptr)
        prof->harvest_ns += step.elapsed_ns();

    assert(
        (!config.stock_on_h || status != TimelineStatus::COMPLETE ||
         harvester.counters.base_descents == base_descents_before) &&
        "invariant 5: a certified component entered the exposed-root-blossom base descent");
    (void)base_descents_before;

    // Back to `G`'s detector ids through the component's own `h_to_det`, which is a strictly
    // ascending subsequence of `H`'s — so a residual sorted inside the component is sorted in `G`,
    // and nothing downstream of harvest has to know the component existed.
    for (uint64_t& defect : harvest.residual) {
        assert(defect < sub.num_nodes());
        defect = sub.h_to_det[defect];
    }
    return status;
}

template <typename HarvestComponent>
Phase1Outcome BallDecoder::decode_impl(
    const std::vector<uint64_t>& dets, BallProfile* prof, bool want_pairs, const HarvestComponent& harvest_component) {
    HiResTimer total_timer;
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

    // ---- §2.2. The connected components of `H`, in ascending root order. This is the *only*
    // partition of the shot: there is no size threshold and no second way to decide a component.
    BallComponentSplit& split = arena.split;
    decompose_ball_components(h, split);

    committed_pair_scratch.clear();
    match_edge_scratch.clear();

    // §M2.9.6 measurement 3. A process-wide counter, so it is read as a delta around the solves.
    uint64_t formations_before = pm::blossom_formation_stats.formations;
    harvester.collect_diagnostics = config.collect_harvest_diagnostics;
    harvester.use_legacy_enumeration = config.use_legacy_harvest_enumeration;

    Phase1Outcome outcome;
    outcome.components_total = (int)split.num_components();

    // ---- §2.2's loop. Every component, whatever its size, is decided by truncated sparse blossom
    // run on it in isolation. A component that truncates escalates the shot on its own; the others
    // are still solved and extracted, because their statuses are what §3's profiler reads and
    // because the escalating shot discards Phase 1 in full anyway (§M3.1).
    HarvestResult part;
    for (size_t c = 0; c < split.num_components(); c++) {
        build_component_subgraph(h, split, c);
        TimelineStatus status = solve_component(prof, harvest_component, part);
        split.status[c] = status == TimelineStatus::TRUNCATED ? BallComponentSplit::COMPONENT_TRUNCATED
                                                              : BallComponentSplit::COMPONENT_COMPLETE;
        if (status == TimelineStatus::TRUNCATED)
            outcome.components_truncated++;
        if (want_pairs)
            drain_component_match_edges();
        accumulate_component_harvest(outcome.harvest, part);
        if (prof != nullptr) {
            // Read off the instance the component was solved on, so it is the largest degree any
            // one solve actually saw. `H`'s own degree distribution is §3.5.2's `degree_hist`.
            for (size_t i = 0; i < arena.split.sub.num_nodes(); i++) {
                prof->max_degree =
                    std::max(prof->max_degree, (int)h_mwpm.mwpm.flooder.graph.nodes[i].neighbors.size());
            }
        }
    }

    // ---- §2.6. **The** escalation predicate, in one place: some component truncated. Everything
    // downstream reads this and nothing re-derives it — not from a residual, not from a certificate,
    // not from a size threshold.
    bool any_component_truncated = outcome.components_truncated > 0;
    outcome.status = any_component_truncated ? TimelineStatus::TRUNCATED : TimelineStatus::COMPLETE;
    assert(any_component_truncated == split.any_truncated() && "§2.6: the status array and the tally disagree");

    // `HarvestResult::residual` is contracted to be sorted ascending in `G`'s ids. Each component's
    // is, but the blocks interleave, so the union is not — this is where that is restored. A no-op
    // on the production path, whose residual is always empty; it is the verification entry points,
    // which harvest every component in full, that actually produce one (§2.5).
    sort_residual(outcome.harvest, harvester.scratch.residual_sort_buffer);

    if (prof != nullptr) {
        prof->intersect_ns = timing.intersect_ns;
        prof->h_build_ns = timing.finalize_ns;
        prof->n_defects = (int)dets.size();
        prof->h_nodes = (int)h.num_nodes();
        prof->h_edges = (int)h.edges.size();
        prof->h_boundary_edges = (int)h.boundary_edges.size();
        // §2.6's tally, and the predicate read off it.
        prof->components_total = outcome.components_total;
        prof->components_truncated = outcome.components_truncated;
        prof->any_component_truncated = any_component_truncated ? 1 : 0;
        // §M7. The shot is certified iff **every** component was: the certificate is a statement
        // about the whole of `H`, and `H`'s matching is the union of the components'.
        if (config.stock_on_h)
            prof->certified = any_component_truncated ? 0 : 1;
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
        }
        prof->shells_materialized = 1;
        prof->restarts = 0;
        prof->blossom_formations = (int)(pm::blossom_formation_stats.formations - formations_before);
        const HarvestResult& result = outcome.harvest;
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
        if (h.num_nodes() != 0)
            prof->mean_degree = 2.0 * (double)h.edges.size() / (double)h.num_nodes();
        prof->total_ns = total_timer.elapsed_ns();
    }
    return outcome;
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

void BallDecoder::verify_decomposition(const Phase1Outcome& outcome, bool full_harvest, bool committed_comparable) {
    // §2.5. The monolithic solve, on a **separate** instance so the two share no state, compared
    // against what the per-component path just produced. This is §1's independence property asked
    // of the machine rather than of the argument for it.
    const BallGraph& h = arena.graph;
    verify_mwpm.rebuild(tables, h, arena, nullptr);

    std::vector<uint64_t>& mono_dets = verify_dets;
    mono_dets.clear();
    mono_dets.reserve(h.num_nodes());
    for (size_t i = 0; i < h.num_nodes(); i++)
        mono_dets.push_back(i);

    TimelineStatus mono_status = TimelineStatus::COMPLETE;
    if (config.stock_on_h) {
        CertificateOutcome certificate = run_stock_and_certify(verify_mwpm.mwpm, mono_dets, horizon, false);
        mono_status = certificate.certified() ? TimelineStatus::COMPLETE : TimelineStatus::TRUNCATED;
    } else {
        mono_status = process_timeline_until_horizon(verify_mwpm.mwpm, mono_dets, horizon);
    }

    auto fail = [](const std::string& what, const std::string& expected, const std::string& actual) {
        throw std::logic_error(
            "The per-component decode disagreed with the monolithic solve on the same H at §2.5 (" + what +
            "): monolithic gave " + expected + ", per-component gave " + actual +
            ". This voids §1's independence property and is a release blocker.");
    };

    // 1. The escalation predicate. This is the one that matters: the set of shots that escalate is a
    //    property of `(DEM, shot, T)` alone (§0), and partitioning `H` must not move it.
    if ((mono_status == TimelineStatus::TRUNCATED) != (outcome.status == TimelineStatus::TRUNCATED)) {
        fail(
            "escalation predicate",
            mono_status == TimelineStatus::TRUNCATED ? "TRUNCATED" : "COMPLETE",
            outcome.status == TimelineStatus::TRUNCATED ? "TRUNCATED" : "COMPLETE");
    }

    HarvestResult mono;
    if (full_harvest) {
        mono = verify_harvester.harvest_to_obs(verify_mwpm.mwpm, mono_dets);
    } else if (mono_status == TimelineStatus::COMPLETE) {
        mono = verify_harvester.extract_only_to_obs(verify_mwpm.mwpm, mono_dets);
    } else {
        // Nothing to compare on a truncated shot the production path never harvested, but the
        // instance still has to be left clean for the next shot.
        abandon_shot(verify_mwpm.mwpm);
        return;
    }

    // 2. On a completing shot: equal weight and equal observable bytes.
    //
    //    `committed` is the **obs** flavour's field. The match-edges flavour appends to the
    //    caller's `CompressedEdge` vector and leaves `committed` at zero on purpose, so on that
    //    flavour this compares `dual_sum_at_truncation` instead — the same flat reduction over the
    //    same live regions, filled by both, and equal to the committed weight on a completed
    //    timeline where every region is frozen and tight.
    //
    //    Above 64 observables the mask is unusable and reads 0 on **both** sides, so the observable
    //    check degenerates to the weight one there. That is stated rather than silently relied on.
    if (outcome.status == TimelineStatus::COMPLETE) {
        if (mono.dual_sum_at_truncation != outcome.harvest.dual_sum_at_truncation)
            fail(
                "dual sum",
                std::to_string(mono.dual_sum_at_truncation),
                std::to_string(outcome.harvest.dual_sum_at_truncation));
        if (committed_comparable) {
            if (mono.committed.weight != outcome.harvest.committed.weight)
                fail(
                    "committed weight",
                    std::to_string(mono.committed.weight),
                    std::to_string(outcome.harvest.committed.weight));
            if (mono.committed.obs_mask != outcome.harvest.committed.obs_mask)
                fail(
                    "committed observables",
                    std::to_string((uint64_t)mono.committed.obs_mask),
                    std::to_string((uint64_t)outcome.harvest.committed.obs_mask));
        }
    }

    // 3. With the full harvest on: the sorted residual set is the union of the per-component
    //    residual sets, and `num_trees` is the sum. The per-component side is already sorted and
    //    already in `G`'s ids; the monolithic side has to be mapped and sorted the same way.
    if (full_harvest) {
        for (uint64_t& defect : mono.residual)
            defect = h.h_to_det[defect];
        std::sort(mono.residual.begin(), mono.residual.end());
        if (mono.residual != outcome.harvest.residual)
            fail("residual set", describe(mono.residual), describe(outcome.harvest.residual));
        if (mono.num_trees != outcome.harvest.num_trees)
            fail("num_trees", std::to_string(mono.num_trees), std::to_string(outcome.harvest.num_trees));
    }
}

HarvestResult BallDecoder::decode_phase1(const std::vector<uint64_t>& dets, BallProfile* prof) {
    if (config.stock_on_h)
        reject_harvest_entry_point();
    // The verification path harvests every component in full, residual and all. It is the same
    // per-component decode the production path takes: §2 leaves exactly one way to decide a
    // component, so there is no longer an "unpruned" variant for an oracle to compare against.
    Phase1Outcome outcome =
        decode_impl(dets, prof, false, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus) {
            return harvester.harvest_to_obs(mwpm, h_dets);
        });
    if (config.verify_component_decomposition)
        verify_decomposition(outcome, true, true);
    if (config.verify_against_g)
        verify_level1(dets, outcome.harvest, nullptr, prof);
    return outcome.harvest;
}

HarvestResult BallDecoder::decode_phase1_to_match_edges(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof) {
    if (config.stock_on_h)
        reject_harvest_entry_point();
    Phase1Outcome outcome =
        decode_impl(dets, prof, true, [this](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus) {
            return harvester.harvest_to_match_edges(mwpm, h_dets, match_edge_scratch);
        });
    committed_pairs = committed_pair_scratch;
    sort_pairs(committed_pairs);

    if (config.verify_component_decomposition)
        verify_decomposition(outcome, true, false);
    if (config.verify_against_g)
        verify_level1(dets, outcome.harvest, &committed_pairs, prof);
    return outcome.harvest;
}

Phase1Outcome BallDecoder::decode_phase1_production(const std::vector<uint64_t>& dets, BallProfile* prof) {
    Phase1Outcome outcome = decode_impl(
        dets, prof, false, [this, prof](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus status) {
            if (status == TimelineStatus::TRUNCATED) {
                // §M3.4's abandon, now per **component** rather than per shot: it is `Mwpm::reset`,
                // which frees the region and node arena pools, so the next component re-allocates
                // them. That is the same teardown an escalating shot always paid, charged once per
                // truncated component instead of once per shot — a shot escalates on ~1 of them, so
                // the multiplier is small, and it is still inside §M3.2's cap on what a cheaper
                // version of this could ever be worth. `abandon_ns` accumulates over the shot,
                // which is why it is `+=`.
                HiResTimer abandon;
                if (prof != nullptr)
                    abandon.start();
                abandon_shot(mwpm);
                if (prof != nullptr)
                    prof->abandon_ns += abandon.elapsed_ns();
                return HarvestResult();
            }
            return harvester.extract_only_to_obs(mwpm, h_dets);
        });
    if (config.verify_component_decomposition)
        verify_decomposition(outcome, false, true);
    return outcome;
}

Phase1Outcome BallDecoder::decode_phase1_production_to_match_edges(
    const std::vector<uint64_t>& dets, std::vector<CommittedPair>& committed_pairs, BallProfile* prof) {
    Phase1Outcome outcome = decode_impl(
        dets, prof, true, [this, prof](pm::Mwpm& mwpm, const std::vector<uint64_t>& h_dets, TimelineStatus status) {
            if (status == TimelineStatus::TRUNCATED) {
                // §M3.4's abandon, now per **component** rather than per shot: it is `Mwpm::reset`,
                // which frees the region and node arena pools, so the next component re-allocates
                // them. That is the same teardown an escalating shot always paid, charged once per
                // truncated component instead of once per shot — a shot escalates on ~1 of them, so
                // the multiplier is small, and it is still inside §M3.2's cap on what a cheaper
                // version of this could ever be worth. `abandon_ns` accumulates over the shot,
                // which is why it is `+=`.
                HiResTimer abandon;
                if (prof != nullptr)
                    abandon.start();
                abandon_shot(mwpm);
                if (prof != nullptr)
                    prof->abandon_ns += abandon.elapsed_ns();
                return HarvestResult();
            }
            return harvester.extract_only_to_match_edges(mwpm, h_dets, match_edge_scratch);
        });
    // §M3.1. The escalating shot re-decodes the whole raw syndrome on `G` and discards **all** of
    // Phase 1 — including the pairs the components that *did* complete contributed, which under §2
    // is most of them. Dropping them here is what keeps a caller from combining them by accident.
    //
    // The verification entry point above deliberately does not do this: it harvests every component
    // in full precisely so that the committed pairs and the residual can be compared against M1 on
    // an escalating shot, which is where §M2.6 level 1 has the most to say.
    if (outcome.status == TimelineStatus::TRUNCATED)
        committed_pair_scratch.clear();
    committed_pairs = committed_pair_scratch;
    sort_pairs(committed_pairs);
    if (config.verify_component_decomposition)
        verify_decomposition(outcome, false, false);
    return outcome;
}

bool BallDecoder::truncated_scheme_escalates(const std::vector<uint64_t>& dets) {
    // The same `H`, built from the same tables at the same `T` filter, and split the same way — the
    // two schemes differ only in the front end, so replaying the decision means replaying the
    // solves, not rebuilding the problem differently.
    compute_seeded_detection_events(dets, seeded_scratch);
    build_ball_graph(tables, seeded_scratch, horizon, arena, config.mode, nullptr, nullptr);
    const BallGraph& h = arena.graph;
    BallComponentSplit& split = arena.split;
    decompose_ball_components(h, split);

    bool escalates = false;
    for (size_t c = 0; c < split.num_components(); c++) {
        build_component_subgraph(h, split, c);
        h_mwpm.rebuild(tables, split.sub, arena, nullptr);

        h_dets_scratch.clear();
        h_dets_scratch.reserve(split.sub.num_nodes());
        for (size_t i = 0; i < split.sub.num_nodes(); i++)
            h_dets_scratch.push_back(i);

        TimelineStatus status = process_timeline_until_horizon(h_mwpm.mwpm, h_dets_scratch, horizon);
        // Nothing is harvested and nothing is read: only the *decision* is wanted. `abandon_shot` is
        // the teardown that works from either outcome, and it is what an escalating shot pays anyway.
        abandon_shot(h_mwpm.mwpm);
        split.status[c] = status == TimelineStatus::TRUNCATED ? BallComponentSplit::COMPONENT_TRUNCATED
                                                              : BallComponentSplit::COMPONENT_COMPLETE;
        escalates = escalates || status == TimelineStatus::TRUNCATED;
    }
    return escalates;
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
