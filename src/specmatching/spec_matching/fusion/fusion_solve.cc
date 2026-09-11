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

#include "specmatching/spec_matching/fusion/fusion_solve.h"

#include <algorithm>
#include <cassert>

#include "specmatching/sparse_blossom/flooder/graph.h"
#include "specmatching/sparse_blossom/flooder/graph_fill_region.h"
#include "specmatching/sparse_blossom/matcher/alternating_tree.h"
#include "specmatching/spec_matching/truncation/harvest.h"

namespace pm {
namespace spec_matching {
namespace fusion {

void EdgeMask::reserve(size_t n_max, size_t entries_max) {
    offsets.reserve(n_max + 1);
    bits.reserve(entries_max);
}

void EdgeMask::rewrite_offsets(const uint32_t* degrees, size_t num_nodes) {
    size_t offsets_capacity = offsets.capacity();
    size_t bits_capacity = bits.capacity();

    offsets.resize(num_nodes + 1);
    uint32_t running = 0;
    for (size_t i = 0; i < num_nodes; i++) {
        offsets[i] = running;
        running += degrees[i];
    }
    offsets[num_nodes] = running;
    // `resize`, not `assign`: the per-piece build writes the mask byte of every entry it writes, so
    // there is nothing here that has to start clear. One byte minimum, so that `bits.data()` is
    // never null on a shot with no edges — a null `edge_mask` is how the flooder is told there is no
    // mask at all, which is emphatically not the same statement.
    bits.resize(std::max<size_t>(running, 1));

    if (offsets.capacity() != offsets_capacity || bits.capacity() != bits_capacity)
        grow_events++;
}

void PieceBuildScratch::reserve(size_t n_max, size_t entries_max) {
    neighbor.reserve(entries_max);
    edge.reserve(entries_max);
    cursor.assign(n_max, 0);
}

void PieceBuildScratch::size_for_shot(size_t entries, size_t num_nodes) {
    neighbor.resize(entries);
    edge.resize(entries);
    if (cursor.size() < num_nodes)
        cursor.resize(num_nodes, 0);
}

void FusionInstance::attach(pm::Mwpm& instance, size_t num_nodes, horizon_int T) {
    mwpm = &instance;
    horizon = T;
    // The offsets are rewritten per shot by `size_mask_for_shot`; all this has to leave behind is a
    // buffer of the right shape and two non-null pointers.
    mask.offsets.assign(num_nodes + 1, 0);
    mask.bits.resize(std::max<size_t>(mask.bits.size(), 1));
    instance.flooder.edge_mask = mask.bits.data();
    instance.flooder.edge_mask_offsets = mask.offsets.data();
    // The cap replaces the global horizon rather than joining it: the clock is shared by every
    // piece of the shot and bounds nothing about any single defect's dual.
    instance.flooder.horizon = pm::NO_HORIZON;
    instance.flooder.dual_cap = T;
    instance.flooder.dual_cap_hit = false;
}

void FusionInstance::size_mask_for_shot(const uint32_t* degrees, size_t num_nodes) {
    assert(mwpm != nullptr && "size_mask_for_shot before attach");
    mask.rewrite_offsets(degrees, num_nodes);
    // Unconditional rather than "only when it grew": a reallocation is exactly the case this is
    // here for, and comparing pointers to decide would be reading a pointer whose buffer may be
    // gone.
    mwpm->flooder.edge_mask = mask.bits.data();
    mwpm->flooder.edge_mask_offsets = mask.offsets.data();
}

void build_piece(
    FusionInstance& fusion,
    PieceBuildScratch& scratch,
    const PieceBuild& in,
    uint32_t piece,
    const uint32_t* defects,
    size_t defect_count,
    const uint32_t* piece_edges,
    size_t piece_edge_count) {
    pm::MatchingGraph& graph = fusion.mwpm->flooder.graph;

    // Step 3 of §3.2, run first: the flooder requires the boundary half-edge to sit at neighbour
    // index 0, so it has to go onto an empty adjacency list. `BallMwpm::rebuild` orders its writes
    // the same way and for the same reason.
    for (size_t at = 0; at < defect_count; at++) {
        uint32_t u = defects[at];
        pm::DetectorNode& node = graph.nodes[u];
        assert(node.neighbors.empty() && "building onto a defect that was never marked unbuilt");
        scratch.cursor[u] = 0;
        if (in.boundary_weight[u] == NO_BOUNDARY_WEIGHT)
            continue;
        node.neighbors.push_back(nullptr);
        node.neighbor_weights.push_back(in.boundary_weight[u]);
        // A dummy carries no observable: it stands for a cut, not for a path to the code boundary.
        node.neighbor_observables.push_back(in.boundary_is_dummy[u] ? (pm::obs_int)0 : in.boundary_obs[u]);
    }

    // Gather. One pass over the piece's edges, each landing in its owner's slice of the scratch. An
    // interior edge writes an entry at both of its endpoints; a crossing edge writes one, at the
    // endpoint that is in this piece — the far endpoint's entry is the far piece's responsibility.
    for (size_t at = 0; at < piece_edge_count; at++) {
        uint32_t e = piece_edges[at];
        const BallGraphEdge& edge = in.edges[e];
        bool interior = in.piece_of[edge.i] == in.piece_of[edge.j];
        if (interior || in.piece_of[edge.i] == piece) {
            uint32_t slot = fusion.mask.offsets[edge.i] + scratch.cursor[edge.i]++;
            scratch.neighbor[slot] = edge.j;
            scratch.edge[slot] = e;
        }
        if (interior || in.piece_of[edge.j] == piece) {
            uint32_t slot = fusion.mask.offsets[edge.j] + scratch.cursor[edge.j]++;
            scratch.neighbor[slot] = edge.i;
            scratch.edge[slot] = e;
        }
    }

    // Order each defect's entries by ascending neighbour, then write them.
    for (size_t at = 0; at < defect_count; at++) {
        uint32_t u = defects[at];
        pm::DetectorNode& node = graph.nodes[u];
        uint32_t base = fusion.mask.offsets[u];
        uint32_t degree = scratch.cursor[u];
        // The boundary half-edge, if there is one, already occupies index 0.
        uint32_t first = (uint32_t)node.neighbors.size();

        // Insertion sort, the design's own suggestion, on `H`'s small degrees. It is the write order
        // and only the write order that this decides: ascending neighbour index is what the
        // adjacency binary searches need, and `(neighbour)` is a strict total order on one defect's
        // entries, so no tie is ever broken.
        for (uint32_t i = 1; i < degree; i++) {
            uint32_t far = scratch.neighbor[base + i];
            uint32_t edge_id = scratch.edge[base + i];
            uint32_t j = i;
            while (j > 0 && scratch.neighbor[base + j - 1] > far) {
                scratch.neighbor[base + j] = scratch.neighbor[base + j - 1];
                scratch.edge[base + j] = scratch.edge[base + j - 1];
                j--;
            }
            scratch.neighbor[base + j] = far;
            scratch.edge[base + j] = edge_id;
        }

        // The boundary half-edge is never masked: it is not a cut, and unmasking is what a fusion
        // does to cuts.
        if (first == 1)
            fusion.mask.set(u, 0, 0);
        for (uint32_t i = 0; i < degree; i++) {
            uint32_t far = scratch.neighbor[base + i];
            uint32_t e = scratch.edge[base + i];
            const BallGraphEdge& edge = in.edges[e];
            node.neighbors.push_back(&graph.nodes[far]);
            node.neighbor_weights.push_back(edge.w_int);
            node.neighbor_observables.push_back(in.edge_obs[e]);

            uint32_t k = first + i;
            if (in.piece_of[far] == piece) {
                fusion.mask.set(u, k, 0);
                continue;
            }
            // A crossing edge is installed masked and stays that way until the fusion at which it
            // becomes interior. That fusion needs to find this entry in O(1), so its index is
            // recorded here, where it is known, rather than binary-searched for later.
            fusion.mask.set(u, k, 1);
            uint32_t r = in.refused_slot[e];
            if (edge.i == u)
                in.cross_index_lo[r] = k;
            else
                in.cross_index_hi[r] = k;
        }
        // The fourth parallel array of the same records; one slot per adjacency entry.
        node.neighbor_implied_weights.resize(node.neighbors.size());
        in.built[u] = 1;
#ifndef NDEBUG
        fusion.debug_node_writes += (uint64_t)node.neighbors.size();
#endif
    }
}

void align_clock_for_new_growth(pm::Mwpm& mwpm) {
    assert(mwpm.flooder.queue.empty() && "the clock can only be nudged between solves, with the queue drained");
    if (mwpm.flooder.queue.cur_time & 1)
        mwpm.flooder.queue.cur_time++;
}

void inject_detection_event(pm::Mwpm& mwpm, horizon_int horizon, uint32_t node_index) {
    (void)horizon;
    pm::GraphFlooder& flooder = mwpm.flooder;
    pm::GraphFillRegion* region = flooder.region_arena.alloc_default_constructed();
    // A default-constructed region is growing with a zero y-intercept, i.e. zero radius at time
    // *zero*. Injected mid-timeline it must have zero radius at the time it is injected.
    region->radius = pm::VaryingCT::growing_varying_with_zero_distance_at_time(flooder.queue.cur_time);
    // `alloc_constructed`, as `Mwpm::create_detection_event` does, so the arena installs the
    // liveness link once the node's lifetime has begun.
    region->alt_tree_node = mwpm.node_arena.alloc_constructed(region);
    flooder.do_region_created_at_empty_detector_node(*region, flooder.graph.nodes[node_index]);
    // A fresh defect sits in no blossom, so its frozen part is zero and its cap is the full `T`.
    flooder.schedule_dual_cap_event(*region, 0);
}

TimelineStatus run_until_settled(pm::Mwpm& mwpm) {
    pm::GraphFlooder& flooder = mwpm.flooder;
    flooder.dual_cap_hit = false;
    while (true) {
        pm::MwpmEvent event = flooder.run_until_next_mwpm_notification();
        if (flooder.dual_cap_hit) {
            // Some defect's dual would have passed `T`. Beyond that point an edge longer than `2T`
            // — one `H` does not contain — could be the next tight edge, so continuing would
            // silently solve a different problem.
            flooder.dual_cap_hit = false;
            return TimelineStatus::TRUNCATED;
        }
        if (event.event_type == pm::NO_EVENT)
            break;
        mwpm.process_event(event);
    }
    // The queue drained. If an alternating tree is still standing, the piece admits no perfect
    // matching at all — odd parity in a piece with no boundary within `T` — which is the same
    // escalation the monolithic truncated path reports, reached by a different route. The test is
    // the O(1) one: tree nodes are checked out of an arena that keeps a liveness vector.
    if (any_alternating_tree_survives(mwpm))
        return TimelineStatus::TRUNCATED;
    // Every region the flooder could reach is matched: to another region, to the real boundary, or
    // to a dummy.
    return TimelineStatus::COMPLETE;
}

TimelineStatus solve_piece(pm::Mwpm& mwpm, horizon_int horizon, const uint32_t* defects, size_t count) {
    align_clock_for_new_growth(mwpm);
    for (size_t i = 0; i < count; i++)
        inject_detection_event(mwpm, horizon, defects[i]);
    return run_until_settled(mwpm);
}

void unmask_crossing_edge(
    FusionInstance& fusion, uint32_t u, uint32_t u_to_v_index, uint32_t v, uint32_t v_to_u_index) {
    fusion.mask.set(u, u_to_v_index, 0);
    fusion.mask.set(v, v_to_u_index, 0);
    pm::GraphFlooder& flooder = fusion.mwpm->flooder;
    flooder.reschedule_events_at_detector_node(flooder.graph.nodes[u]);
    flooder.reschedule_events_at_detector_node(flooder.graph.nodes[v]);
}

bool set_boundary_and_release(pm::Mwpm& mwpm, uint32_t u, pm::weight_int weight, pm::obs_int obs_mask) {
    pm::GraphFlooder& flooder = mwpm.flooder;
    pm::DetectorNode& node = flooder.graph.nodes[u];
    assert(!node.neighbors.empty() && node.neighbors[0] == nullptr && "defect has no boundary half-edge to raise");
    node.neighbor_weights[0] = weight;
    node.neighbor_observables[0] = obs_mask;

    pm::GraphFillRegion* region = node.region_that_arrived_top;
    // Matched to the boundary *through this defect*: `match.region == nullptr` is the boundary
    // match, and `match.edge.loc_from` names the defect whose boundary edge it went through. A
    // region matched to the boundary through some other defect keeps that match; the dummy that is
    // going away is not the one holding it.
    bool matched_here = region != nullptr && region->alt_tree_node == nullptr && region->match.region == nullptr &&
                        region->match.edge.loc_from == &node;
    if (!matched_here) {
        // The weight changed under a node that stays where it is; its next-event calculation has to
        // be redone. (Released regions get this for free from `set_region_growing`.)
        flooder.reschedule_events_at_detector_node(node);
        return false;
    }

    // Release. Detaching the match and making the region an alternating-tree root is exactly what
    // `handle_tree_hitting_match` does to the far side of a match it steals, and `set_region_growing`
    // resumes the radius rather than resetting it — `then_growing_at_time` keeps the distance it has
    // right now and only changes the slope.
    region->match.clear();
    region->alt_tree_node = mwpm.node_arena.alloc_constructed(region);
    flooder.set_region_growing(*region);
    return true;
}

pm::MatchingResult extract_root(pm::Mwpm& mwpm, const uint32_t* defects, size_t count) {
    pm::MatchingResult result;
    for (size_t i = 0; i < count; i++) {
        pm::DetectorNode& node = mwpm.flooder.graph.nodes[defects[i]];
        pm::GraphFillRegion* region = node.region_that_arrived_top;
        // Already extracted: the shatter that took this defect's region reset every node it owned.
        if (region == nullptr)
            continue;
        assert(region->alt_tree_node == nullptr && "extracting a root that still has an alternating tree");
        assert(
            (region->match.region != nullptr || region->match.edge.loc_from != nullptr) &&
            "a settled region is neither matched nor in a tree");
        result += mwpm.shatter_blossom_and_extract_matches(region);
    }
    return result;
}

bool feasibility_holds(
    const FusionInstance& fusion, const uint32_t* defects, size_t count, const char** failure_out) {
    const pm::Mwpm& mwpm = *fusion.mwpm;
    const pm::MatchingGraph& graph = mwpm.flooder.graph;
    pm::cumulative_time_int now = mwpm.flooder.queue.cur_time;

    auto fail = [&](const char* why) {
        if (failure_out != nullptr)
            *failure_out = why;
        return false;
    };

    for (size_t i = 0; i < count; i++) {
        const pm::DetectorNode& node = graph.nodes[defects[i]];
        pm::cumulative_time_int dual = node.local_radius().get_distance_at_time(now);
        if (dual > fusion.horizon)
            return fail("a defect's dual exceeds T");
        if (dual < 0)
            return fail("a defect's dual is negative");
    }

    for (size_t i = 0; i < count; i++) {
        uint32_t index = defects[i];
        const pm::DetectorNode& node = graph.nodes[index];
        pm::cumulative_time_int rad1 = node.local_radius().get_distance_at_time(now);
        for (size_t k = 0; k < node.neighbors.size(); k++) {
            if (fusion.mask.get(index, (uint32_t)k))
                continue;  // Masked: not an edge of the problem currently being solved.
            const pm::DetectorNode* neighbor = node.neighbors[k];
            pm::cumulative_time_int weight = node.neighbor_weights[k];
            if (neighbor == nullptr) {
                if (rad1 > weight)
                    return fail("a boundary edge has negative slack");
                continue;
            }
            // Two nodes of one region have grown into each other by construction; the flooder's own
            // collision test skips them for the same reason.
            if (node.has_same_owner_as(*neighbor))
                continue;
            pm::cumulative_time_int rad2 = neighbor->local_radius().get_distance_at_time(now);
            if (rad1 + rad2 > weight)
                return fail("an unmasked edge has negative slack");
        }
    }
    return true;
}

}  // namespace fusion
}  // namespace spec_matching
}  // namespace pm
