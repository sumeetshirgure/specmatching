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

#ifndef PYREMATCHING_TWO_PHASE_TRUNCATION_HARVEST_H
#define PYREMATCHING_TWO_PHASE_TRUNCATION_HARVEST_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "pyrematching/two_phase/truncation/exposed_blossom.h"
#include "pyrematching/two_phase/truncation/horizon.h"

namespace pm {
namespace two_phase {

/// One committed alternating-tree pair, recorded for the debug-build tight-edge invariant:
/// `d(inner, outer) == Y(inner) + Y(outer)`. Filled only when diagnostics are requested; the
/// distance itself is not computed here (it needs a search-graph Dijkstra, which the tests run).
struct TightPairRecord {
    uint64_t inner_node;
    uint64_t outer_node;
    pm::total_weight_int inner_dual_sum;
    pm::total_weight_int outer_dual_sum;
};

/// What one shot's truncated timeline yielded.
///
/// The commit policy is exact and deliberate: **the residual is exactly the exposed defects, one
/// per surviving alternating tree.** Everything else — every frozen match, and every
/// (inner_region, outer_region) pair of every non-root tree node — is committed. A policy that
/// dumped whole trees into the residual would break the separation property the two-phase design
/// rests on.
struct HarvestResult {
    /// The exposed defects, sorted ascending. `residual.size() == num_trees`.
    std::vector<uint64_t> residual;
    /// `Y(u)` for each residual defect, aligned with `residual`. Must equal the horizon exactly;
    /// this is the separation invariant, and a violation voids the Phase-2 error bound.
    std::vector<pm::total_weight_int> residual_dual_sum;

    /// The committed matching, obs flavour. Unused by the match-edges flavour, which appends to
    /// the caller's `match_edges` vector instead.
    pm::MatchingResult committed;

    /// Pairs already matched when the horizon was reached.
    int committed_pairs_frozen{0};
    /// Pairs committed by walking the surviving trees (non-root `inner`/`outer` pairs).
    int committed_pairs_tree{0};
    /// Pairs committed by pairing the even cycle around the base of an exposed root blossom.
    int committed_pairs_blossom_cycle{0};
    /// Regions committed to the boundary.
    int committed_boundary{0};

    int num_trees{0};
    /// Largest surviving tree, counted in `AltTreeNode`s (each of which is an inner/outer pair,
    /// except the root, which is just the exposed outer region).
    ///
    /// **Zero unless `Harvester::collect_diagnostics` is set.** §M2.9.1's direct enumeration
    /// classifies each node from its own fields and never learns which tree a node belongs to, so
    /// this now costs a root id per node — a parent chase, which is the thing that enumeration
    /// exists to remove. It is kept because §M2.9.6 measurement 1 needs it, but it is paid for only
    /// when asked for.
    int largest_tree_size{0};
    /// How many surviving tree roots were blossoms rather than plain regions. These are the
    /// highest-risk case: the residual defect has to be the blossom's base.
    int exposed_root_blossoms{0};

    /// `Sum_S y_S` over every region alive at truncation. A lower bound on the exact MWPM optimum,
    /// and the denominator of the per-shot certificate ratio.
    pm::total_weight_int dual_sum_at_truncation{0};
    /// `max_S y_S` over the same regions. Compared against the cluster weight-diameter in the M1
    /// exit artifact: it is the quantity that says whether a single dual can outgrow the horizon.
    pm::total_weight_int max_region_dual{0};

    /// §M2.9.6 measurements. All of these are filled only when `Harvester::collect_diagnostics` is
    /// set, and are reported in counts and depths rather than nanoseconds — the point is a
    /// cycle-count model that does not depend on this machine's cache hierarchy.
    ///
    /// Deepest blossom nesting alive at truncation: 0 when no blossom survived, 1 for a blossom of
    /// plain regions, 2 once a blossom contains a blossom. Measurement 2, and the quantity §M2.9.4's
    /// choice of base-descent mechanism turns on.
    int max_blossom_nesting_depth{0};
    /// Most members (leaf regions, i.e. detection events) in any one surviving blossom.
    int max_blossom_members{0};
    /// The same two, restricted to *exposed root* blossoms — the only ones the base descent of
    /// §M2.9.4 has to walk. These are what decide that section, and they are not the same
    /// distribution as the matched blossoms above.
    int max_exposed_blossom_depth{0};
    int max_exposed_blossom_members{0};
    /// Top-level matched regions that are blossoms, and so have to be shattered to yield a pairing.
    /// Both endpoints of a blossom-to-blossom match count. Measurement for §M2.9.5.
    int matched_blossom_shatters{0};
    /// Modelled serial depth of harvest, in dependent operations. See `harvest.cc` for the model.
    int harvest_dependent_depth{0};

    /// Wall time of harvest's four stages, for the software A/B of §M2.9.6 ("how much of harvest's
    /// 7–17% was enumeration rather than extraction"). Filled only under `collect_diagnostics`, and
    /// deliberately not part of the H1 identity comparison.
    long long enumerate_ns{0};
    long long reduce_ns{0};
    long long base_descent_ns{0};
    long long shatter_ns{0};

    void clear();

    /// Everything H1 compares: the output proper plus the counters, but not the timings and not
    /// `harvest_dependent_depth` (a model, not an output).
    bool identical_output_to(const HarvestResult& other) const;
    std::string describe_difference(const HarvestResult& other) const;
};

/// Reusable scratch buffers, so that harvesting a shot allocates nothing in steady state.
struct HarvestScratch {
    /// One representative per committed match: the endpoint that extraction is driven from.
    std::vector<pm::GraphFillRegion*> top_regions;
    /// A copy of the node arena's live vector, taken before the sweep that deletes from it.
    std::vector<pm::AltTreeNode*> live_nodes;
    std::vector<pm::AltTreeNode*> tree_roots;
    std::vector<std::pair<uint64_t, pm::total_weight_int>> residual_sort_buffer;

    /// Used by the M1.3 oracle path only (`Harvester::use_legacy_enumeration`), and by the debug
    /// equivalence check of §M2.9.2.
    std::unordered_set<pm::GraphFillRegion*> seen_regions;
    std::unordered_set<pm::AltTreeNode*> seen_roots;
    std::vector<pm::AltTreeNode*> tree_walk_stack;
    /// Diagnostics only: tree sizes keyed by root, for `largest_tree_size`.
    std::unordered_map<pm::AltTreeNode*, int> tree_sizes;

    void clear();
};

/// What harvest actually did, counted rather than timed (§M3.3 X9).
///
/// The §M3.4 bypass is a claim about which *code* runs, not about how long it takes, so it is
/// checked as one: on the production path `tree_nodes_visited` and `base_descents` must be
/// identically zero, on every shot, and no amount of scheduler noise can make that pass by
/// accident. Cumulative across shots; the caller resets when it wants a per-shot reading.
struct HarvestCounters {
    /// Live regions the enumeration walked. This is extraction's own trip count, not overhead.
    uint64_t regions_visited{0};
    /// Live `AltTreeNode`s the enumeration walked. Zero on the extract-only path by construction.
    uint64_t tree_nodes_visited{0};
    /// Calls into M1.4's base descent — the design's highest-risk item, and the one §M3.4 exists to
    /// take off the production path.
    uint64_t base_descents{0};
    /// Matched top-level regions handed to the shatter/extract step.
    uint64_t extractions{0};
    /// Calls to the full harvest and to the extract-only path respectively, so a test can tell
    /// which one a decoder actually took.
    uint64_t full_harvests{0};
    uint64_t extract_only_harvests{0};

    void reset() {
        *this = HarvestCounters();
    }
};

/// Extracts the committed matching and the residual from a truncated (or completed) timeline.
///
/// Owns the scratch buffers, so keep one alive across shots. `tight_pairs_out`, when non-null,
/// collects a `TightPairRecord` per tree-committed pair; it is a diagnostic and costs nothing when
/// left null.
struct Harvester {
    HarvestScratch scratch;
    std::vector<TightPairRecord>* tight_pairs_out{nullptr};
    HarvestCounters counters;

    /// Runs the M1.3 enumeration — sweep the shot's detection events, dedup regions with a visited
    /// stamp, find each tree's root and descend `children` from it — instead of §M2.9.1's direct
    /// enumeration. Kept compiled in as the H1 oracle, exactly as `verify_against_g` is for M2.
    /// The commit *policy* is identical either way; only the enumeration differs.
    bool use_legacy_enumeration{false};

    /// Fills the §M2.9.6 measurements, including `largest_tree_size`, which direct enumeration no
    /// longer gets for free (a root id per node is a parent chase, which is exactly what §M2.9.1
    /// removes). Off by default so the hot path pays nothing; the M1/M2/M2.9 artifacts and the
    /// identity tests turn it on.
    bool collect_diagnostics{false};

    /// Harvest into an observable mask plus a weight. Requires at most 64 observables.
    HarvestResult harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

    /// Harvest into compressed match edges, appended to `match_edges`. Works for any number of
    /// observables; the caller turns the edges into observables with
    /// `Mwpm::extract_paths_from_match_edges`.
    HarvestResult harvest_to_match_edges(
        pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges);

    /// §M3.4's production path: **extraction only**.
    ///
    /// Valid only on a `COMPLETE` timeline, where no alternating tree survives. That is exactly the
    /// case in which the two jobs this drops — classifying `AltTreeNode`s into tree commits, and
    /// M1.4's base descent through an exposed root blossom — have nothing to do; a `TRUNCATED`
    /// timeline escalates instead (§M3.1) and its Phase-1 result is discarded, so the jobs would be
    /// computing something that is never read.
    ///
    /// The result carries the committed observables and weight, the frozen/boundary commit counts
    /// and `dual_sum_at_truncation` — the last because the §M4.2 certificate reads it and it is a
    /// flat reduction over the same live-region walk that extraction already does. `residual`,
    /// `num_trees`, `largest_tree_size` and `exposed_root_blossoms` are zero, and are zero as
    /// *facts* about a completed timeline rather than as unfilled fields.
    ///
    /// **The full harvest above stays compiled in and is the oracle** (§M3.3 X8, and §M2.6 level 1
    /// still runs against it): deleting it would delete what licenses this.
    HarvestResult extract_only_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);
    HarvestResult extract_only_to_match_edges(
        pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges);

   private:
    template <typename ExtractMatched>
    HarvestResult extract_only_impl(pm::Mwpm& mwpm, const ExtractMatched& extract_matched);
    template <typename ExtractMatched, typename ExtractExposed>
    HarvestResult harvest_impl(
        pm::Mwpm& mwpm,
        const std::vector<uint64_t>& detection_events,
        const ExtractMatched& extract_matched,
        const ExtractExposed& extract_exposed);
    template <typename ExtractMatched, typename ExtractExposed>
    HarvestResult harvest_impl_legacy(
        pm::Mwpm& mwpm,
        const std::vector<uint64_t>& detection_events,
        const ExtractMatched& extract_matched,
        const ExtractExposed& extract_exposed);
};

/// Throws the shot away without harvesting it: §M3.4's `TRUNCATED` branch.
///
/// An escalating shot re-decodes the whole syndrome on `G` and discards Phase 1's partial result
/// (§M3.1), so nothing that harvest would extract is ever read — but the instance still has to be
/// left clean for the next shot, and the surviving regions and tree nodes are what harvest would
/// otherwise have handed back to their arenas.
///
/// This is `Mwpm::reset` — the same teardown the stock driver runs when it cannot find a perfect
/// matching — and it is the expensive kind: it sweeps `H`'s node vector and frees the arena pools,
/// so the next shot re-allocates them. That is deliberate. Writing a cheaper "return everything
/// without extracting it" would duplicate the extraction's ownership rules, which is the one part
/// of this code with a silent failure mode; and the cost is paid on `q <= 3.25e-4` of shots, over
/// `H`'s node count rather than `G`'s. Both bounds are measured, and §M3.2 caps what any cleverer
/// version of this could ever be worth at ~0.02% amortised.
inline void abandon_shot(pm::Mwpm& mwpm) {
    mwpm.reset();
}

/// Does any alternating tree survive? **O(1)** — one `empty()` on the node arena's live vector,
/// which §M2.9.1 maintains inside `Arena<T>` itself.
///
/// This is the branch §M3.4's production control flow turns on, and the reason it has to be O(1) is
/// that an O(n) sweep over regions to decide whether to skip an O(n) harvest defeats the point.
/// `any_alternating_tree_survives_by_sweep` is the same question answered the slow way, and the two
/// are asserted equal on every shot in debug builds (§M3.3 X9).
inline bool any_alternating_tree_survives(const pm::Mwpm& mwpm) {
    return !mwpm.node_arena.live.empty();
}

/// The full sweep the O(1) branch above is checked against: a region is in a tree iff it carries an
/// `alt_tree_node`. Tests and debug asserts only.
bool any_alternating_tree_survives_by_sweep(const pm::Mwpm& mwpm);

/// Puts the instance back into the state the next shot's timeline preamble expects. Called at the
/// end of every harvest; exposed for tests and for callers that drive the pieces themselves.
///
/// This is `Mwpm::reset` minus the arena teardown, which the post-harvest state does not need and
/// which would otherwise dominate the per-shot cost.
void reset_for_next_shot(pm::Mwpm& mwpm);

/// Convenience wrappers that allocate their own scratch. Prefer a long-lived `Harvester` on any
/// path that runs more than one shot.
HarvestResult harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);
HarvestResult harvest_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_TRUNCATION_HARVEST_H
