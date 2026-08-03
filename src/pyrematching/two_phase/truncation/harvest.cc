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

#include "pyrematching/two_phase/truncation/harvest.h"

#include <algorithm>
#include <cassert>

#include "pyrematching/sparse_blossom/flooder/graph.h"
#include "pyrematching/sparse_blossom/matcher/alternating_tree.h"
#include "pyrematching/sparse_blossom/tracker/flood_check_event.h"

namespace pm {
namespace two_phase {

void HarvestResult::clear() {
    residual.clear();
    residual_dual_sum.clear();
    committed = pm::MatchingResult();
    committed_pairs_frozen = 0;
    committed_pairs_tree = 0;
    committed_pairs_blossom_cycle = 0;
    committed_boundary = 0;
    num_trees = 0;
    largest_tree_size = 0;
    exposed_root_blossoms = 0;
    dual_sum_at_truncation = 0;
    max_region_dual = 0;
}

void HarvestScratch::clear() {
    top_regions.clear();
    seen_regions.clear();
    tree_roots.clear();
    seen_roots.clear();
    tree_walk_stack.clear();
    residual_sort_buffer.clear();
}

template <typename ExtractMatched, typename ExtractExposed>
HarvestResult Harvester::harvest_impl(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    const ExtractMatched& extract_matched,
    const ExtractExposed& extract_exposed) {
    HarvestResult result;
    scratch.clear();

    auto& flooder = mwpm.flooder;
    auto& nodes = flooder.graph.nodes;
    const pm::DetectorNode* node_base = nodes.data();
    // `process_timeline_until_horizon` has already advanced the clock to exactly the horizon, so
    // every dual read below is the dual at `T`.
    const pm::cumulative_time_int time = flooder.queue.cur_time;

    // The negative-weight detection events are extracted alongside the shot's own, exactly as the
    // stock decode path does: they are real detection events whose regions have to be resolved.
    const std::vector<uint64_t>* event_lists[2] = {&detection_events, &flooder.negative_weight_detection_events};

    // ---- Phase A: classify. Reads only; nothing is mutated until every region is accounted for.
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            pm::GraphFillRegion* top = node.region_that_arrived_top;
            if (!scratch.seen_regions.insert(top).second)
                continue;
            scratch.top_regions.push_back(top);
        }
    }

    int frozen_half_pairs = 0;
    for (pm::GraphFillRegion* top : scratch.top_regions) {
        // Every live region is `top` or nested inside it, so this covers the whole dual solution.
        top->do_op_for_each_descendant_and_self([&](pm::GraphFillRegion* region) {
            pm::total_weight_int dual = region->radius.get_distance_at_time(time);
            result.dual_sum_at_truncation += dual;
            result.max_region_dual = std::max(result.max_region_dual, dual);
        });

        if (top->alt_tree_node != nullptr) {
            auto* root = const_cast<pm::AltTreeNode*>(top->alt_tree_node->find_root());
            if (scratch.seen_roots.insert(root).second)
                scratch.tree_roots.push_back(root);
        } else if (top->match.region != nullptr) {
            frozen_half_pairs++;
        } else {
            // Not in a tree and not matched to a region: matched to the boundary. An unmatched
            // region outside every tree would mean the primal state is inconsistent.
            assert(top->match.edge.loc_from != nullptr && "top-level region is neither matched nor in a tree");
            result.committed_boundary++;
        }
    }
    assert(frozen_half_pairs % 2 == 0 && "a frozen match should contribute both of its endpoints");
    result.committed_pairs_frozen = frozen_half_pairs / 2;

    // ---- Phase B: commit the surviving trees and expose one defect per tree.
    for (pm::AltTreeNode* root : scratch.tree_roots) {
        result.num_trees++;

        int tree_size = 0;
        scratch.tree_walk_stack.clear();
        scratch.tree_walk_stack.push_back(root);
        while (!scratch.tree_walk_stack.empty()) {
            pm::AltTreeNode* tree_node = scratch.tree_walk_stack.back();
            scratch.tree_walk_stack.pop_back();
            tree_size++;
            if (tree_node->inner_region != nullptr) {
                // Committed, tight: the inner/outer pair of a non-root node is already a matched
                // pair of the current primal, joined by a tight edge.
                result.committed_pairs_tree++;
                if (tight_pairs_out != nullptr) {
                    const pm::CompressedEdge& edge = tree_node->inner_to_outer_edge;
                    tight_pairs_out->push_back(
                        TightPairRecord{
                            (uint64_t)(edge.loc_from - node_base),
                            (uint64_t)(edge.loc_to - node_base),
                            nested_dual_sum(*edge.loc_from, time),
                            nested_dual_sum(*edge.loc_to, time)});
                }
            }
            for (const pm::AltTreeEdge& child : tree_node->children)
                scratch.tree_walk_stack.push_back(child.alt_tree_node);
        }
        result.largest_tree_size = std::max(result.largest_tree_size, tree_size);

        pm::GraphFillRegion* exposed = root->outer_region;
        if (!exposed->blossom_children.empty())
            result.exposed_root_blossoms++;

        // Matches and freezes every non-root inner/outer pair, and dismantles the tree. The root
        // is left alone beyond clearing its back-pointer: its outer region is the residual.
        mwpm.shatter_descendants_into_matches_and_freeze(*root);

        pm::DetectorNode* base_node = find_exposed_base_node(*exposed, time);
        assert(base_node != nullptr);
        pm::total_weight_int base_dual_sum = nested_dual_sum(*base_node, time);

        uint64_t exposed_defect = 0;
        result.committed_pairs_blossom_cycle += (int)extract_exposed(*exposed, exposed_defect);
        assert(exposed_defect == (uint64_t)(base_node - node_base));

        result.residual.push_back(exposed_defect);
        result.residual_dual_sum.push_back(base_dual_sum);
    }

    // ---- Phase C: extract everything that is matched. Shattering resets the detector nodes it
    // consumes, so re-reading `region_that_arrived` each time is what stops a pair being extracted
    // twice from its two endpoints. This mirrors the stock extraction helpers.
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            extract_matched(node.region_that_arrived_top);
        }
    }

    // The residual is emitted sorted ascending, with its dual sums kept aligned. Sorted through a
    // scratch buffer so the only allocations per shot are the result's own.
    scratch.residual_sort_buffer.clear();
    scratch.residual_sort_buffer.reserve(result.residual.size());
    for (size_t i = 0; i < result.residual.size(); i++)
        scratch.residual_sort_buffer.emplace_back(result.residual[i], result.residual_dual_sum[i]);
    std::sort(scratch.residual_sort_buffer.begin(), scratch.residual_sort_buffer.end());
    for (size_t i = 0; i < scratch.residual_sort_buffer.size(); i++) {
        result.residual[i] = scratch.residual_sort_buffer[i].first;
        result.residual_dual_sum[i] = scratch.residual_sort_buffer[i].second;
    }
    assert(result.residual.size() == (size_t)result.num_trees);

    reset_for_next_shot(mwpm);
    return result;
}

void reset_for_next_shot(pm::Mwpm& mwpm) {
    // Every region and alternating tree node has already been handed back to its arena by the
    // extraction above, so — unlike `Mwpm::reset` — the arenas are left alone. Destroying them
    // would free the pools and make the next shot malloc every region again, which costs more than
    // the whole truncated timeline.
    assert(
        mwpm.node_arena.allocated.size() == mwpm.node_arena.available.size() &&
        "harvest left an alternating tree node checked out");
    assert(
        mwpm.flooder.region_arena.allocated.size() == mwpm.flooder.region_arena.available.size() &&
        "harvest left a graph fill region checked out");

    // Shattering already reset every node the shot owned, via `cleanup_shell_area`, and a node a
    // shrinking region released was zeroed on the way out. What a truncated timeline leaves behind
    // that a completed one does not is the queue: it drains itself when the matching completes,
    // but here it can still hold exempt shrink events. A node with an event still in the queue is
    // exactly a node whose tracker still believes it has one queued, and that belief would stop
    // the next shot from ever enqueueing for that node.
    //
    // So clear those trackers and drop the queue — work proportional to what is left over, rather
    // than to the size of the graph. `Mwpm::reset`'s whole-graph sweep would cost more per shot
    // than the timeline it is cleaning up after.
    for (auto& bucket : mwpm.flooder.queue.bit_buckets) {
        for (auto& event : bucket) {
            if (event.tentative_event_type == pm::LOOK_AT_NODE)
                event.data_look_at_node->node_event_tracker.clear();
        }
    }
    mwpm.flooder.queue.clear();

    // Shrink-event trackers live in graph fill regions, which are default-constructed afresh when
    // the arena hands them out again, so they need no cleaning here.

#ifndef NDEBUG
    // The reasoning above is subtle enough to be worth checking outright in debug builds: the
    // instance must be indistinguishable from a fresh one.
    for (const auto& node : mwpm.flooder.graph.nodes) {
        assert(node.region_that_arrived == nullptr);
        assert(node.region_that_arrived_top == nullptr);
        assert(node.reached_from_source == nullptr);
        assert(node.radius_of_arrival == 0);
        assert(node.wrapped_radius_cached == 0);
        assert(node.observables_crossed_from_source == 0);
        assert(!node.node_event_tracker.has_queued_time);
    }
#endif
}

HarvestResult Harvester::harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    HarvestResult result;
    pm::MatchingResult committed;
    result = harvest_impl(
        mwpm,
        detection_events,
        [&](pm::GraphFillRegion* region) {
            committed += mwpm.shatter_blossom_and_extract_matches(region);
        },
        [&](pm::GraphFillRegion& exposed, uint64_t& exposed_defect_out) {
            return shatter_exposed_blossom_and_extract_matches(
                mwpm, exposed, mwpm.flooder.queue.cur_time, exposed_defect_out, committed);
        });
    result.committed = committed;
    return result;
}

HarvestResult Harvester::harvest_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges) {
    return harvest_impl(
        mwpm,
        detection_events,
        [&](pm::GraphFillRegion* region) {
            mwpm.shatter_blossom_and_extract_match_edges(region, match_edges);
        },
        [&](pm::GraphFillRegion& exposed, uint64_t& exposed_defect_out) {
            return shatter_exposed_blossom_and_extract_match_edges(
                mwpm, exposed, mwpm.flooder.queue.cur_time, exposed_defect_out, match_edges);
        });
}

HarvestResult harvest_to_obs(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    Harvester harvester;
    return harvester.harvest_to_obs(mwpm, detection_events);
}

HarvestResult harvest_to_match_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<pm::CompressedEdge>& match_edges) {
    Harvester harvester;
    return harvester.harvest_to_match_edges(mwpm, detection_events, match_edges);
}

}  // namespace two_phase
}  // namespace pm
