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

    void clear();
};

/// Reusable scratch buffers, so that harvesting a shot allocates nothing in steady state.
struct HarvestScratch {
    std::vector<pm::GraphFillRegion*> top_regions;
    std::unordered_set<pm::GraphFillRegion*> seen_regions;
    std::vector<pm::AltTreeNode*> tree_roots;
    std::unordered_set<pm::AltTreeNode*> seen_roots;
    std::vector<pm::AltTreeNode*> tree_walk_stack;
    std::vector<std::pair<uint64_t, pm::total_weight_int>> residual_sort_buffer;

    void clear();
};

/// Extracts the committed matching and the residual from a truncated (or completed) timeline.
///
/// Owns the scratch buffers, so keep one alive across shots. `tight_pairs_out`, when non-null,
/// collects a `TightPairRecord` per tree-committed pair; it is a diagnostic and costs nothing when
/// left null.
struct Harvester {
    HarvestScratch scratch;
    std::vector<TightPairRecord>* tight_pairs_out{nullptr};

    /// Harvest into an observable mask plus a weight. Requires at most 64 observables.
    HarvestResult harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

    /// Harvest into compressed match edges, appended to `match_edges`. Works for any number of
    /// observables; the caller turns the edges into observables with
    /// `Mwpm::extract_paths_from_match_edges`.
    HarvestResult harvest_to_match_edges(
        pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges);

   private:
    template <typename ExtractMatched, typename ExtractExposed>
    HarvestResult harvest_impl(
        pm::Mwpm& mwpm,
        const std::vector<uint64_t>& detection_events,
        const ExtractMatched& extract_matched,
        const ExtractExposed& extract_exposed);
};

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
