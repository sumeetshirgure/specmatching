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

#ifndef SPECMATCHING_SPEC_MATCHING_TRUNCATION_EXPOSED_BLOSSOM_H
#define SPECMATCHING_SPEC_MATCHING_TRUNCATION_EXPOSED_BLOSSOM_H

#include <vector>

#include "specmatching/sparse_blossom/matcher/mwpm.h"
#include "specmatching/spec_matching/truncation/horizon.h"

namespace pm {
namespace spec_matching {

/// `Y(u)`: the sum of `y_S` over every region `S` containing the detection event `u`, evaluated at
/// `time`.
///
/// This is the dual value that the separation invariant constrains: a residual defect harvested at
/// horizon `T` must have `Y(u) == T` exactly, because an exposed region grows at rate 1 from time
/// 0 and is never frozen. `Y(u) <= T` holds for every detection event, since at any moment exactly
/// one region containing `u` is varying and it varies at rate at most 1.
///
/// Computed from `compute_wrapped_radius()` rather than the `wrapped_radius_cached` field so that
/// it does not depend on the cache being in sync.
inline pm::cumulative_time_int nested_dual_sum(const pm::DetectorNode& node, pm::cumulative_time_int time) {
    if (node.region_that_arrived_top == nullptr)
        return 0;
    return (node.region_that_arrived_top->radius + node.compute_wrapped_radius()).get_distance_at_time(time);
}

/// Finds the detection event that a region is *exposed at*: the source node whose nested dual sum
/// is maximal, recursed down through however many levels of nested blossoms there are.
///
/// For a region that is not a blossom this is simply the detection event the region grew from. For
/// a blossom it is the base of the odd cycle — the member that the alternating path enters
/// through, whose dual has grown continuously and therefore equals the current time.
///
/// The maximum is not always unique. On a degenerate graph a region can be matched and pulled back
/// into a tree at the same instant, spending zero time frozen, so its dual also tracks the clock.
/// Any such member is a legitimate answer: every edge of the odd cycle is tight, so removing any
/// member leaves a tight perfect matching on the rest, and the separation invariant asks only for
/// `Y(u) == T`. Ties break towards the lower detector index, for determinism.
///
/// Does not mutate anything: callers must resolve the base *before* shattering starts, because
/// shattering detaches sub-blossoms and destroys the nesting the answer is derived from.
pm::DetectorNode* find_exposed_base_node(pm::GraphFillRegion& region, pm::cumulative_time_int time);

/// Shatters an *exposed* (unmatched) region — the outer region of a surviving alternating tree
/// root — into the one defect that stays exposed plus a perfect matching on everything else.
///
/// `pm::Mwpm::shatter_blossom_and_extract_matches` and friends all assume a *matched* blossom:
/// they read `region->match.edge.loc_from` to decide which sub-blossom is the base. An exposed
/// blossom has no match edge, so the base is instead identified by `find_exposed_base_node`, and
/// the remaining even cycle around it is paired and extracted exactly as
/// `pair_and_shatter_subblossoms_and_extract_matches` does.
///
/// Unlike the matched case, the region's own `y_S` is *not* added to the extracted weight: the dual
/// of an exposed region belongs to the residual defect it is exposed at, not to any committed pair.
///
/// Consumes the region (it and all its sub-blossoms are returned to the arena, and their shell
/// areas are reset), writes the index of the exposed detection event to `exposed_defect_out`, and
/// returns the number of pairs it committed.
size_t shatter_exposed_blossom_and_extract_matches(
    pm::Mwpm& mwpm,
    pm::GraphFillRegion& root_blossom,
    pm::cumulative_time_int time,
    uint64_t& exposed_defect_out,
    pm::MatchingResult& res);

/// Match-edges flavour of `shatter_exposed_blossom_and_extract_matches`, for graphs with more than
/// 64 observables (where the `obs_int` mask is not usable).
size_t shatter_exposed_blossom_and_extract_match_edges(
    pm::Mwpm& mwpm,
    pm::GraphFillRegion& root_blossom,
    pm::cumulative_time_int time,
    uint64_t& exposed_defect_out,
    std::vector<pm::CompressedEdge>& match_edges);

}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_TRUNCATION_EXPOSED_BLOSSOM_H
