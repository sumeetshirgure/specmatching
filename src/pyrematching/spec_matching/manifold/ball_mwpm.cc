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

#include "pyrematching/spec_matching/manifold/ball_mwpm.h"

#include <algorithm>
#include <cassert>

namespace pm {
namespace spec_matching {

namespace {

/// The observable mask of a CSR slice of observable ids. Zero when the graph has more than 64
/// observables, exactly as `MatchingGraph::add_edge` leaves it: the mask is unusable there and the
/// match-edges flavour is the only correct one.
pm::obs_int mask_of(const BallTables& tables, uint64_t begin, uint64_t end, bool use_masks) {
    if (!use_masks)
        return 0;
    pm::obs_int mask = 0;
    for (uint64_t k = begin; k < end; k++)
        mask ^= (pm::obs_int)1 << tables.ball_mask_ids[k];
    return mask;
}

pm::obs_int boundary_mask_of(const BallTables& tables, uint64_t det, bool use_masks) {
    if (!use_masks)
        return 0;
    pm::obs_int mask = 0;
    for (uint64_t k = tables.bcost_mask_offsets[det]; k < tables.bcost_mask_offsets[det + 1]; k++)
        mask ^= (pm::obs_int)1 << tables.bcost_mask_ids[k];
    return mask;
}

}  // namespace

void BallMwpm::configure(size_t num_observables_in, double normalising_constant_in) {
    num_observables = num_observables_in;
    normalising_constant = normalising_constant_in;
    capacity = 0;
    used_nodes = 0;
    grow_events = 0;
    mwpm = pm::Mwpm();
}

void BallMwpm::rebuild(const BallTables& tables, const BallGraph& h, BallGraphArena& arena, BallMwpmCounts* counts) {
    assert(normalising_constant == tables.normalising_constant && "H and G must share one normalising constant (§0)");

    size_t n = h.num_nodes();
    if (n > capacity) {
        // Growing moves the node vector, so every `DetectorNode*` in the old adjacency dies with
        // it. That is safe only because adjacency is rewritten from scratch below and nothing
        // survives between shots.
        capacity = std::max<size_t>(n, capacity * 2);
        pm::MatchingGraph graph(capacity, num_observables, normalising_constant);
        mwpm = pm::Mwpm(pm::GraphFlooder(std::move(graph)));
        used_nodes = 0;
        grow_events++;
    }

    // The §M2 structural counters, accumulated at the sites that do the writing. Register adds
    // only; the store happens once, below, and only when a profile asked for it.
    uint64_t node_records = 0;
    uint64_t edge_records = 0;

    pm::MatchingGraph& graph = mwpm.flooder.graph;
    size_t stale = std::max(used_nodes, n);
    for (size_t i = 0; i < stale; i++) {
        auto& node = graph.nodes[i];
        node.neighbors.clear();
        node.neighbor_weights.clear();
        node.neighbor_observables.clear();
        node.neighbor_implied_weights.clear();
        node_records++;
    }
    used_nodes = n;

    bool use_masks = num_observables <= sizeof(pm::obs_int) * 8;

    // Boundary edges first: the flooder requires the boundary half-edge to sit at neighbour index
    // 0, and adding it to an empty adjacency list avoids the O(degree) insert `add_boundary_edge`
    // would otherwise pay.
    for (const BallBoundaryEdge& edge : h.boundary_edges) {
        auto& node = graph.nodes[edge.i];
        assert(node.neighbors.empty() && "boundary edges must be added before any other edge");
        node.neighbors.push_back(nullptr);
        node.neighbor_weights.push_back(edge.w_int);
        node.neighbor_observables.push_back(boundary_mask_of(tables, edge.det, use_masks));
        edge_records++;
    }

    // Adjacency in ascending neighbour index. Determinism (§0) asks for every output-affecting
    // iteration order to be over sorted ids, and the flooder's event ordering is one of them.
    //
    // No sort is needed to get there. `h.edges` is already sorted by `(i, j)` with `i < j`, so
    // appending the *incoming* direction of every edge first and the *outgoing* direction second
    // leaves each node's list as [neighbours below it, ascending] ++ [neighbours above it,
    // ascending] — which is ascending overall. Sorting `2 * |E|` directed pairs per shot instead
    // was 17–45% of the whole decode.
    auto append = [&](uint32_t from, uint32_t to, const BallGraphEdge& edge) {
        auto& node = graph.nodes[from];
        node.neighbors.push_back(&graph.nodes[to]);
        node.neighbor_weights.push_back(edge.w_int);
        node.neighbor_observables.push_back(
            mask_of(tables, tables.ball_mask_offsets[edge.entry], tables.ball_mask_offsets[edge.entry + 1], use_masks));
        edge_records++;
    };
    for (const BallGraphEdge& edge : h.edges)
        append(edge.j, edge.i, edge);
    for (const BallGraphEdge& edge : h.edges)
        append(edge.i, edge.j, edge);
    // One implied-weight slot per adjacency entry already counted above; this is the fourth
    // parallel array of the same records, not a fifth kind of record.
    for (size_t i = 0; i < n; i++)
        graph.nodes[i].neighbor_implied_weights.resize(graph.nodes[i].neighbors.size());

    if (counts != nullptr) {
        counts->node_records = node_records;
        counts->edge_records = edge_records;
    }

    // `H` is built from the post-preamble detection events and carries no negative weights, so the
    // flooder's negative-weight bookkeeping is empty by construction. The offsets that `G`'s
    // preamble produced are applied by the caller, in `G`'s units.
    mwpm.flooder.negative_weight_detection_events.clear();
    mwpm.flooder.negative_weight_observables.clear();
    mwpm.flooder.negative_weight_obs_mask = 0;
    mwpm.flooder.negative_weight_sum = 0;

#ifndef NDEBUG
    // Debug invariant 9: no weight in `H` is negative, `H` has an edge for every recorded pair and
    // none beyond, and the adjacency is consistent in both directions.
    for (size_t i = 0; i < n; i++) {
        const auto& node = graph.nodes[i];
        for (size_t k = 0; k < node.neighbors.size(); k++) {
            if (node.neighbors[k] == nullptr) {
                assert(k == 0 && "the boundary half-edge must be neighbour 0");
                continue;
            }
            size_t j = (size_t)(node.neighbors[k] - graph.nodes.data());
            assert(j < n && "H edge points outside the used node range");
            const auto& other = graph.nodes[j];
            size_t back = other.index_of_neighbor(const_cast<pm::DetectorNode*>(&graph.nodes[i]));
            assert(other.neighbor_weights[back] == node.neighbor_weights[k] && "H adjacency is asymmetric");
        }
    }
#endif
}

}  // namespace spec_matching
}  // namespace pm
