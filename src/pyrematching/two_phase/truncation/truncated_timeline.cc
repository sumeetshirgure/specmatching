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

#include "pyrematching/two_phase/truncation/truncated_timeline.h"

#include <algorithm>
#include <cassert>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/sparse_blossom/flooder/graph.h"
#include "pyrematching/sparse_blossom/flooder/graph_fill_region.h"
#include "pyrematching/two_phase/truncation/harvest.h"

namespace pm {
namespace two_phase {

void TimelineDepthModel::begin(const pm::MatchingGraph& graph) {
    node_depth.assign(graph.nodes.size(), 0);
    node_base = graph.nodes.data();
    depth = 0;
    events = 0;
}

void TimelineDepthModel::observe(const pm::MwpmEvent& event) {
    if (node_base == nullptr || event.event_type == pm::NO_EVENT)
        return;
    events++;

    auto index_of = [&](const pm::DetectorNode* node) -> size_t {
        if (node == nullptr)
            return SIZE_MAX;
        size_t index = (size_t)(node - node_base);
        return index < node_depth.size() ? index : SIZE_MAX;
    };

    // Collect the nodes this event depends on, take the deepest, and write one deeper back to all
    // of them. `touched` is a small stack buffer for the common case; a blossom shatter falls back
    // to a sweep of the blossom's area.
    size_t touched[2] = {SIZE_MAX, SIZE_MAX};
    int deepest = 0;

    switch (event.event_type) {
        case pm::REGION_HIT_REGION: {
            const pm::CompressedEdge& edge = event.region_hit_region_event_data.edge;
            touched[0] = index_of(edge.loc_from);
            touched[1] = index_of(edge.loc_to);
            break;
        }
        case pm::REGION_HIT_BOUNDARY: {
            const pm::CompressedEdge& edge = event.region_hit_boundary_event_data.edge;
            touched[0] = index_of(edge.loc_from);
            touched[1] = index_of(edge.loc_to);
            break;
        }
        case pm::BLOSSOM_SHATTER: {
            // A shatter names regions rather than an edge, and it un-does a merge, so it depends on
            // everything the blossom had absorbed. Sweeping its area is the faithful reading, and
            // this path only runs under measurement.
            pm::GraphFillRegion* blossom = event.blossom_shatter_event_data.blossom_region;
            std::vector<size_t> members;
            blossom->do_op_for_each_node_in_total_area([&](pm::DetectorNode* node) {
                size_t index = index_of(node);
                if (index != SIZE_MAX) {
                    members.push_back(index);
                    deepest = std::max(deepest, node_depth[index]);
                }
            });
            for (size_t index : members)
                node_depth[index] = std::max(node_depth[index], deepest + 1);
            depth = std::max(depth, deepest + 1);
            return;
        }
        default:
            return;
    }

    for (size_t index : touched) {
        if (index != SIZE_MAX)
            deepest = std::max(deepest, node_depth[index]);
    }
    for (size_t index : touched) {
        if (index != SIZE_MAX)
            node_depth[index] = deepest + 1;
    }
    depth = std::max(depth, deepest + 1);
}

namespace {

/// The timeline proper. `depth_model` is null on the hot path, and the only difference the measured
/// path makes is one virtual-free call per event.
TimelineStatus run_timeline(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    horizon_int horizon,
    TimelineDepthModel* depth_model) {
    // Reuses the stock preamble verbatim: empty-queue precondition, cur_time = 0, and the
    // negative-weight detection event marking dance.
    pm::begin_timeline(mwpm, detection_events);

    // Exception-safe: the flooder goes back to the sentinel horizon however we leave this scope,
    // so the stock decode path never observes a finite horizon.
    HorizonGuard guard(mwpm.flooder, horizon);

    while (true) {
        auto event = mwpm.flooder.run_until_next_mwpm_notification();
        if (event.event_type == pm::NO_EVENT)
            break;
        if (depth_model != nullptr)
            depth_model->observe(event);
        mwpm.process_event(event);
    }

    // Alternating tree nodes are allocated from `node_arena` and returned to it when a tree is
    // resolved, so "some node is still checked out" is exactly "some alternating tree survives".
    //
    // §M3.4 requires this decision to be O(1), because it is the branch that decides whether to run
    // the harvest at all and an O(n) sweep to skip an O(n) pass buys nothing. It is: §M2.9.1's
    // arena liveness vector answers it with one `empty()`. Debug builds cross-check the same
    // question asked two slower ways — the free-list comparison this used to read, and a full sweep
    // over live regions looking for an `alt_tree_node` (§M3.3 X9).
    bool trees_survive = any_alternating_tree_survives(mwpm);
    assert(trees_survive == (mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size()));
    assert(trees_survive == any_alternating_tree_survives_by_sweep(mwpm));
    if (!trees_survive)
        return TimelineStatus::COMPLETE;

    if (horizon != pm::NO_HORIZON) {
        // Advance the clock to exactly `T`. The queue holds nothing at or before the horizon
        // (node events past it were never enqueued, and the loop above stopped at the first
        // exempt shrink event past it), so the surviving regions grow undisturbed to `T`. Doing
        // this here, once, is what makes `Y(u) == T` an exact equality at harvest.
        mwpm.flooder.queue.cur_time = horizon;
    }
    return TimelineStatus::TRUNCATED;
}

}  // namespace

TimelineStatus process_timeline_until_horizon(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, horizon_int horizon) {
    return run_timeline(mwpm, detection_events, horizon, nullptr);
}

TimelineStatus process_timeline_until_horizon_measured(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    horizon_int horizon,
    TimelineDepthModel& depth_model) {
    depth_model.begin(mwpm.flooder.graph);
    return run_timeline(mwpm, detection_events, horizon, &depth_model);
}

}  // namespace two_phase
}  // namespace pm
