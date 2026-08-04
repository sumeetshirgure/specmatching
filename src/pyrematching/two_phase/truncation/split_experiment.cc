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

#include "pyrematching/two_phase/truncation/split_experiment.h"

#include <algorithm>
#include <cassert>

#include "pyrematching/sparse_blossom/flooder/graph.h"
#include "pyrematching/sparse_blossom/flooder/graph_fill_region.h"

namespace pm {
namespace two_phase {

namespace {

constexpr size_t NOT_A_CHILD = SIZE_MAX;

/// Index, among `blossom.blossom_children`, of the child whose nesting contains `node`.
///
/// Upstream reads this as `region->match.edge.loc_from->region_that_arrived_top`, which works there
/// because by that point the children have had `clear_blossom_parent_ignoring_wrapped_radius()`
/// called and each one *is* the top for its own nodes. Nothing is detached here, so the same answer
/// is obtained by climbing the `blossom_parent` chain until the parent is the blossom itself.
size_t find_child_index(const pm::GraphFillRegion& blossom, const pm::DetectorNode* node) {
    if (node == nullptr || node->region_that_arrived == nullptr)
        return NOT_A_CHILD;
    const pm::GraphFillRegion* current = node->region_that_arrived;
    while (current != nullptr && current->blossom_parent != &blossom)
        current = current->blossom_parent;
    if (current == nullptr)
        return NOT_A_CHILD;
    for (size_t i = 0; i < blossom.blossom_children.size(); i++) {
        if (blossom.blossom_children[i].region == current)
            return i;
    }
    return NOT_A_CHILD;
}

size_t choose_split(const pm::GraphFillRegion& blossom, const pm::DetectorNode* attach, SplitRule rule) {
    size_t index = find_child_index(blossom, attach);
    // Under `ROTATED` the inherited match edge no longer attaches inside the sub-blossom that
    // inherited it, one level down, so the attach point genuinely is outside this cycle. Index 0 is
    // then as good a starting point as any; the rotation below still moves it.
    if (index == NOT_A_CHILD)
        index = 0;
    if (rule == SplitRule::ROTATED)
        index = (index + 1) % blossom.blossom_children.size();
    return index;
}

/// Counts the blossoms in a region's nesting whose tight odd cycle has a non-zero total observable
/// mask. A cycle that is homologically trivial flips every observable an even number of times, so
/// its edges XOR to zero; this is that statement, counted rather than argued.
void count_cycle_obs(const pm::GraphFillRegion& region, uint64_t& examined, uint64_t& nonzero) {
    if (region.blossom_children.empty())
        return;
    pm::obs_int cycle = 0;
    for (const pm::RegionEdge& child : region.blossom_children) {
        cycle ^= child.edge.obs_mask;
        count_cycle_obs(*child.region, examined, nonzero);
    }
    examined++;
    if (cycle != 0)
        nonzero++;
}

int nesting_depth(const pm::GraphFillRegion& region) {
    if (region.blossom_children.empty())
        return 0;
    int deepest = 0;
    for (const pm::RegionEdge& child : region.blossom_children)
        deepest = std::max(deepest, nesting_depth(*child.region));
    return deepest + 1;
}

/// `y_S` of a region, exactly the quantity upstream's extraction accumulates.
inline pm::total_weight_int y_of(const pm::GraphFillRegion& region) {
    return const_cast<pm::GraphFillRegion&>(region).radius.y_intercept();
}

}  // namespace

SplitEvaluation evaluate_shatter(
    const pm::GraphFillRegion& region,
    const pm::GraphFillRegion* partner,
    const pm::CompressedEdge& edge,
    SplitRule rule) {
    // The two base cases of `shatter_blossom_and_extract_matches`, in the same order.
    if (partner != nullptr && region.blossom_children.empty() && partner->blossom_children.empty())
        return SplitEvaluation{edge.obs_mask, y_of(region) + y_of(*partner)};
    if (partner == nullptr && region.blossom_children.empty())
        return SplitEvaluation{edge.obs_mask, y_of(region)};

    SplitEvaluation result;
    const pm::GraphFillRegion* inner_region = &region;
    const pm::GraphFillRegion* inner_partner = partner;

    // `pair_and_shatter_subblossoms_and_extract_matches`, on this side of the match: the blossom's
    // own dual is committed, the cycle is split, and the even path that is left pairs up along it.
    // The sub-blossom at the split point inherits the external match and is handled by the recursive
    // call at the bottom.
    if (!region.blossom_children.empty()) {
        result.weight += y_of(region);
        const auto& children = region.blossom_children;
        size_t count = children.size();
        size_t split = choose_split(region, edge.loc_from, rule);
        for (size_t i = 0; i + 1 < count; i += 2) {
            const pm::RegionEdge& first = children[(split + i + 1) % count];
            const pm::RegionEdge& second = children[(split + i + 2) % count];
            SplitEvaluation sub = evaluate_shatter(*first.region, second.region, first.edge, rule);
            result.obs_mask ^= sub.obs_mask;
            result.weight += sub.weight;
        }
        inner_region = children[split].region;
    }

    // And the same on the partner's side. Upstream splits it at `region->match.edge`'s far end,
    // which is `edge.loc_to`.
    if (inner_partner != nullptr && !inner_partner->blossom_children.empty()) {
        result.weight += y_of(*inner_partner);
        const auto& children = inner_partner->blossom_children;
        size_t count = children.size();
        size_t split = choose_split(*inner_partner, edge.loc_to, rule);
        for (size_t i = 0; i + 1 < count; i += 2) {
            const pm::RegionEdge& first = children[(split + i + 1) % count];
            const pm::RegionEdge& second = children[(split + i + 2) % count];
            SplitEvaluation sub = evaluate_shatter(*first.region, second.region, first.edge, rule);
            result.obs_mask ^= sub.obs_mask;
            result.weight += sub.weight;
        }
        inner_partner = children[split].region;
    }

    SplitEvaluation rest = evaluate_shatter(*inner_region, inner_partner, edge, rule);
    result.obs_mask ^= rest.obs_mask;
    result.weight += rest.weight;
    return result;
}

SplitEvaluation evaluate_shatter(const pm::GraphFillRegion& region, SplitRule rule) {
    return evaluate_shatter(region, region.match.region, region.match.edge, rule);
}

void probe_split_dependence(pm::Mwpm& mwpm, SplitExperimentStats& stats) {
    stats.shots++;
    bool saw_matched_blossom = false;

    for (pm::GraphFillRegion* region : mwpm.flooder.region_arena.live) {
        if (region->blossom_parent != nullptr || region->alt_tree_node != nullptr)
            continue;  // Not top-level, or still in a surviving tree.
        pm::GraphFillRegion* partner = region->match.region;
        if (partner == nullptr)
            continue;  // Boundary match: no cycle, nothing to split.
        // One probe per pair, from the endpoint holding the lower detector id.
        if (!(region->match.edge.loc_from < region->match.edge.loc_to))
            continue;
        if (region->blossom_children.empty() && partner->blossom_children.empty())
            continue;  // No blossom on either side, so the split point does not exist.

        saw_matched_blossom = true;
        stats.splittable_matches++;
        stats.max_nesting_depth =
            std::max(stats.max_nesting_depth, std::max(nesting_depth(*region), nesting_depth(*partner)));

        SplitEvaluation upstream = evaluate_shatter(*region, SplitRule::UPSTREAM);
        SplitEvaluation rotated = evaluate_shatter(*region, SplitRule::ROTATED);
        if (upstream.obs_mask != rotated.obs_mask)
            stats.obs_differed++;
        if (upstream.weight != rotated.weight)
            stats.weight_differed++;
        // What "do not shatter this blossom at all" would report: the region-level match edge, and
        // nothing from the interior.
        if (upstream.obs_mask != region->match.edge.obs_mask)
            stats.interior_obs_nonzero++;
        count_cycle_obs(*region, stats.cycles_examined, stats.cycles_with_nonzero_obs);
        count_cycle_obs(*partner, stats.cycles_examined, stats.cycles_with_nonzero_obs);
    }

    if (saw_matched_blossom)
        stats.shots_with_matched_blossoms++;
}

}  // namespace two_phase
}  // namespace pm
