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

#ifndef PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_MWPM_H
#define PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_MWPM_H

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "pyrematching/spec_matching/manifold/ball_graph.h"
#include "pyrematching/spec_matching/manifold/ball_tables.h"

namespace pm {
namespace spec_matching {

/// What a rebuild of the `pm::Mwpm` on `H` actually writes, counted at the write sites (§M2
/// structural counters). Filled only when a non-null pointer is passed.
///
/// `node_records` is the number of `DetectorNode`s whose adjacency state this rebuild resets — not
/// `H`'s node count: the reset has to cover the *previous* shot's nodes too, so it is
/// `max(used_nodes, |H|)`.
///
/// `edge_records` is the number of directed adjacency entries appended: adjacency is stored
/// directed, so an undirected edge of `H` is two of them, and a boundary edge is one. Each record
/// is four parallel-array element writes (`neighbors`, `neighbor_weights`, `neighbor_observables`,
/// `neighbor_implied_weights`), the last of which is written by the resize pass over the nodes.
struct BallMwpmCounts {
    uint64_t node_records{0};
    uint64_t edge_records{0};

    inline uint64_t total() const {
        return node_records + edge_records;
    }
};

/// A `pm::Mwpm` living on the per-shot graph `H`, rebuilt in place from shot to shot.
///
/// **The discretisation trap (§M2.4).** `UserGraph::to_mwpm` recomputes a normalising constant from
/// its own weight range. Applied to `H` — whose weights are already integer `G`-distances — it would
/// rescale everything and destroy the correspondence with both `G`'s metric and `T`. So `H`'s
/// matching graph is built directly from the pre-discretised integers, and it inherits `G`'s
/// `normalising_constant` verbatim rather than deriving one.
///
/// The node pool is grown, never shrunk, and adjacency is rewritten in place, so a steady-state
/// rebuild allocates nothing.
struct BallMwpm {
    pm::Mwpm mwpm;

    /// How many nodes the underlying matching graph holds; `>= h.num_nodes()` after any rebuild.
    size_t capacity{0};
    /// How many of them the last rebuild used, so the next one only has to clear those.
    size_t used_nodes{0};
    size_t num_observables{0};
    double normalising_constant{0};
    /// Rebuilds that had to enlarge the node pool. Stops increasing once a corpus's largest shot
    /// has been seen; that is the steady state debug invariant 11 asks for.
    uint64_t grow_events{0};

    /// Points the instance at a detector graph's observable count and normalising constant. Must be
    /// called before the first `rebuild`.
    void configure(size_t num_observables, double normalising_constant);

    /// Rewrites `mwpm`'s matching graph to be `h`. `arena` supplies the scratch for canonicalising
    /// the adjacency, so that this allocates nothing in steady state.
    void rebuild(const BallTables& tables, const BallGraph& h, BallGraphArena& arena, BallMwpmCounts* counts = nullptr);
};

}  // namespace spec_matching
}  // namespace pm

#endif  // PYREMATCHING_SPEC_MATCHING_MANIFOLD_BALL_MWPM_H
