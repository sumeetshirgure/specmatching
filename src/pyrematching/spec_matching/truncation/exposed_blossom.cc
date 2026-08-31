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

#include "pyrematching/spec_matching/truncation/exposed_blossom.h"

#include <algorithm>
#include <cassert>
#include <limits>

#include "pyrematching/sparse_blossom/flooder/graph.h"

namespace pm {
namespace spec_matching {

pm::DetectorNode* find_exposed_base_node(pm::GraphFillRegion& region, pm::cumulative_time_int time) {
    pm::DetectorNode* best = nullptr;
    pm::cumulative_time_int best_dual_sum = std::numeric_limits<pm::cumulative_time_int>::min();

    region.do_op_for_each_descendant_and_self([&](pm::GraphFillRegion* r) {
        // Only the leaves of the blossom nesting own a detection event; the blossoms above them
        // are aggregates. A leaf's source is the first node it ever owned, and shell areas only
        // ever shrink from the back, so it is still at the front.
        if (!r->blossom_children.empty() || r->shell_area.empty())
            return;
        pm::DetectorNode* source = r->shell_area[0];
        pm::cumulative_time_int dual_sum = nested_dual_sum(*source, time);
        if (best == nullptr || dual_sum > best_dual_sum || (dual_sum == best_dual_sum && source < best)) {
            best = source;
            best_dual_sum = dual_sum;
        }
    });

    return best;
}

namespace {

/// Shared body of the two flavours. `extract_pair` is handed the first region of each matched pair
/// once the pair has been formed, and is responsible for extracting and freeing it.
template <typename ExtractPair>
size_t shatter_exposed_blossom_impl(
    pm::Mwpm& mwpm,
    pm::GraphFillRegion& root_blossom,
    pm::cumulative_time_int time,
    uint64_t& exposed_defect_out,
    const ExtractPair& extract_pair) {
    pm::DetectorNode* base_node = find_exposed_base_node(root_blossom, time);
    assert(base_node != nullptr && "an exposed region always owns at least one detection event");
    exposed_defect_out = (uint64_t)(base_node - &mwpm.flooder.graph.nodes[0]);

    // The chain of regions from the exposed root down to the leaf that owns the base defect. It is
    // resolved up front because shattering detaches sub-blossoms level by level, destroying the
    // very nesting the answer is read from.
    std::vector<pm::GraphFillRegion*> chain;
    for (pm::GraphFillRegion* r = base_node->region_that_arrived; r != nullptr; r = r->blossom_parent)
        chain.push_back(r);
    std::reverse(chain.begin(), chain.end());
    assert(!chain.empty() && chain.front() == &root_blossom);

    size_t pairs_committed = 0;
    for (size_t level = 0; level + 1 < chain.size(); level++) {
        pm::GraphFillRegion* blossom = chain[level];
        pm::GraphFillRegion* base_child = chain[level + 1];

        // Nodes the blossom itself reached after it formed belong to no surviving region.
        blossom->cleanup_shell_area();
        // Hand each sub-blossom back its own nodes, so that the ordinary (matched) shatter
        // routines can recurse into whichever of them are blossoms themselves.
        for (auto& child : blossom->blossom_children)
            child.region->clear_blossom_parent_ignoring_wrapped_radius();

        size_t num_children = blossom->blossom_children.size();
        size_t base_index = 0;
        for (size_t i = 0; i < num_children; i++) {
            if (blossom->blossom_children[i].region == base_child) {
                base_index = i;
                break;
            }
        }

        // Removing the base from the odd cycle leaves an even path, which pairs up along the cycle
        // exactly as in `pair_and_shatter_subblossoms_and_extract_matches`.
        for (size_t i = 0; i + 1 < num_children; i += 2) {
            auto& re1 = blossom->blossom_children[(base_index + i + 1) % num_children];
            auto& re2 = blossom->blossom_children[(base_index + i + 2) % num_children];
            re1.region->add_match(re2.region, re1.edge);
            extract_pair(re1.region);
            pairs_committed++;
        }

        // Deliberately no `res.weight += blossom->radius.y_intercept()` here, unlike the matched
        // case: an exposed region's dual is part of `Y(base)` and belongs to the residual defect,
        // not to any pair committed in this shot.
        mwpm.flooder.region_arena.del(blossom);
    }

    pm::GraphFillRegion* base_leaf = chain.back();
    base_leaf->cleanup_shell_area();
    mwpm.flooder.region_arena.del(base_leaf);
    return pairs_committed;
}

}  // namespace

size_t shatter_exposed_blossom_and_extract_matches(
    pm::Mwpm& mwpm,
    pm::GraphFillRegion& root_blossom,
    pm::cumulative_time_int time,
    uint64_t& exposed_defect_out,
    pm::MatchingResult& res) {
    return shatter_exposed_blossom_impl(mwpm, root_blossom, time, exposed_defect_out, [&](pm::GraphFillRegion* region) {
        res += mwpm.shatter_blossom_and_extract_matches(region);
    });
}

size_t shatter_exposed_blossom_and_extract_match_edges(
    pm::Mwpm& mwpm,
    pm::GraphFillRegion& root_blossom,
    pm::cumulative_time_int time,
    uint64_t& exposed_defect_out,
    std::vector<pm::CompressedEdge>& match_edges) {
    return shatter_exposed_blossom_impl(mwpm, root_blossom, time, exposed_defect_out, [&](pm::GraphFillRegion* region) {
        mwpm.shatter_blossom_and_extract_match_edges(region, match_edges);
    });
}

}  // namespace spec_matching
}  // namespace pm
