// Copyright 2022 PyMatching Contributors
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

#ifndef PYREMATCHING2_GRAPH_FILL_REGION_H
#define PYREMATCHING2_GRAPH_FILL_REGION_H

#include <vector>

#include "pyrematching/sparse_blossom/arena.h"
#include "pyrematching/sparse_blossom/flooder/match.h"
#include "pyrematching/sparse_blossom/flooder_matcher_interop/region_edge.h"
#include "pyrematching/sparse_blossom/flooder_matcher_interop/varying.h"
#include "pyrematching/sparse_blossom/tracker/queued_event_tracker.h"

namespace pm {

class AltTreeNode;

struct GraphFillRegion {
    /// If this region has merged with others into a blossom, this is that blossom.
    GraphFillRegion* blossom_parent;
    /// The topmost fill region that contains this region. This field must be kept up to date as
    /// the region structure changes.
    GraphFillRegion* blossom_parent_top;
    /// If this is a top-level region (not a blossom child), this is the alternating tree that
    /// it is part of. Note that it may be a degenerate alternating tree with just a single
    /// graph fill region (this one).
    pm::AltTreeNode* alt_tree_node;
    /// How much this region has grown since its creation (and whether it is currently growing).
    /// For graph fill regions starting from detection events, this is just the actual radius
    /// of the region. For graph fill regions starting from blossom-creation events, this is how
    /// much the blossom has grown since it was created (it's the *extra* radius, not the total
    /// radius starting from the detection events).
    pm::VaryingCT radius;
    /// Event tracker for shrink events. Handles ensuring at least, and ideally exactly, one event
    /// to look at the region is in the event queue.
    QueuedEventTracker shrink_event_tracker;

    /// Slot in `Arena<GraphFillRegion>::live`, maintained by the arena (pyrematching §M2.9.2).
    /// Declaring it is what opts this type into arena liveness tracking; nothing in the flooder or
    /// the matcher reads or writes it. It lets harvest enumerate every region alive at truncation —
    /// and hence, by one local `blossom_parent == nullptr` test, every *top-level* region — without
    /// sweeping the shot's detection events and deduplicating them against a visited stamp.
    ///
    /// **Declared here rather than at the end of the struct on purpose.** `shrink_event_tracker`
    /// ends at offset 44 and `match` is 8-aligned, so there is a four-byte hole here that this
    /// field fits in exactly; appended at the end it would instead grow `GraphFillRegion` from 128
    /// bytes to 136, and on this machine that showed up as a measurable regression in the stock
    /// decode path. Keep it here.
    ///
    /// The move constructor deliberately leaves it at its default: a region that is moved out of
    /// the arena is a new object that is not in the arena's live vector, and the arena never moves
    /// the regions it owns.
    uint32_t arena_live_index = pm::ARENA_NOT_LIVE;

    /// If the region is matched, as opposed to growing/shrinking, this says what it is matched to.
    pm::Match match;

    /// If this region is a blossom, these are its child regions along with the cyclic paths
    /// between the children.
    std::vector<pm::RegionEdge> blossom_children;
    /// The set of nodes directly owned by this region. Note that this vector does not include
    /// the nodes indirectly owned by this region that are owned by the blossom children of this
    /// region (or their children or etc).
    std::vector<pm::DetectorNode*> shell_area;

    void cleanup_shell_area();

    GraphFillRegion();
    GraphFillRegion(GraphFillRegion&&);
    GraphFillRegion(const GraphFillRegion&) = delete;
    bool tree_equal(const pm::GraphFillRegion& other) const;

    void add_match(pm::GraphFillRegion* match, const pm::CompressedEdge& edge);

    template <typename Callable>
    void do_op_for_each_node_in_total_area(const Callable& func);
    template <typename Callable>
    void do_op_for_each_descendant_and_self(const Callable& func);
    void clear_blossom_parent();
    void clear_blossom_parent_ignoring_wrapped_radius();
    void wrap_into_blossom(GraphFillRegion* new_blossom_parent_and_top);

    /// Determines if rhs is an ancestor of, or the same as, lhs.
    /// Caution: The equality used here is reference equality (not tree equality like operator==).
    bool operator<=(const GraphFillRegion& rhs) const;
    /// Determines if rhs is an ancestor of, but not the same as, lhs.
    bool operator<(const GraphFillRegion& rhs) const;
    /// Determines if rhs is a descendent of, or the same as, lhs.
    /// Caution: The equality used here is reference equality (not tree equality like operator==).
    bool operator>=(const GraphFillRegion& rhs) const;
    /// Determines if rhs is a descendent of, but not the same as, lhs.
    bool operator>(const GraphFillRegion& rhs) const;

    bool operator==(const GraphFillRegion& rhs) const;
    bool operator!=(const GraphFillRegion& rhs) const;
};

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_node_in_total_area(const Callable& func) {
    for (size_t i = 0; i < shell_area.size(); i++) {
        func(shell_area[shell_area.size() - i - 1]);
    }
    for (auto& child : blossom_children) {
        child.region->do_op_for_each_node_in_total_area(func);
    }
}

template <typename Callable>
inline void pm::GraphFillRegion::do_op_for_each_descendant_and_self(const Callable& func) {
    func(this);
    for (RegionEdge& child : blossom_children) {
        child.region->do_op_for_each_descendant_and_self(func);
    }
}

}  // namespace pm

#endif  // PYREMATCHING2_GRAPH_FILL_REGION_H
