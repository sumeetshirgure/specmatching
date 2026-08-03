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

#ifndef PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_GRAPH_H
#define PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_GRAPH_H

#include <cstdint>
#include <vector>

#include "pyrematching/two_phase/manifold/ball_tables.h"

namespace pm {
namespace two_phase {

/// How the shot-time ball intersection is performed. Both are built, benchmarked and kept (§M2.3);
/// they must produce byte-identical graphs (test B5).
enum class BallGraphBuildMode {
    /// Walk `v`'s ball entries and test the syndrome bit of each. Cost `n * |B|`.
    SCAN,
    /// AND `v`'s ball bitset against the syndrome bitset window and extract the set bits. Cost
    /// `n * ball_word_len[v]` words.
    BITSET,
};

struct BallGraphEdge {
    uint32_t i;
    uint32_t j;
    pm::weight_int w_int;
    /// Index into the ball pools of the canonical path behind this edge, always looked up from the
    /// *lower* detector id of the pair. The canonical tie-break is not symmetric, so fixing the
    /// direction here is what makes the edge's observable set a function of the pair alone.
    uint64_t entry;
};

struct BallBoundaryEdge {
    uint32_t i;
    pm::weight_int w_int;
    /// Detector id of the node, i.e. the index into the `bcost_*` tables.
    uint64_t det;
};

/// The per-shot graph `H` of §M2.1: nodes are the shot's (post-preamble) detection events, indexed
/// `0..n-1` in ascending detector id; edges are every defect pair within `2T`, weighted by the
/// **exact** integer distance in `G`'s discretised metric; boundary edges are every defect with
/// `bcost <= T`. Never persisted.
struct BallGraph {
    /// Detector id of each `H` node, ascending.
    std::vector<uint64_t> h_to_det;
    /// Sorted by `(i, j)` with `i < j`, so `H` does not depend on which build mode produced it.
    std::vector<BallGraphEdge> edges;
    /// Sorted by `i`.
    std::vector<BallBoundaryEdge> boundary_edges;

    inline size_t num_nodes() const {
        return h_to_det.size();
    }

    void clear() {
        h_to_det.clear();
        edges.clear();
        boundary_edges.clear();
    }
};

/// Owns everything `H`'s construction needs across shots, so that steady-state building allocates
/// nothing (debug invariant 11). `det_to_h` and the syndrome bitset are sized once to the detector
/// graph and cleared through touched lists, never by re-zeroing.
struct BallGraphArena {
    BallGraph graph;

    std::vector<uint32_t> det_to_h;
    std::vector<uint64_t> syndrome_words;
    std::vector<uint32_t> touched_words;

    /// Counts builds that had to grow a buffer. After a warmup shot on a representative corpus this
    /// must stop increasing; that is what "zero allocations in steady state" is asserted as.
    uint64_t grow_events{0};

    void reset_for_graph(size_t num_detector_nodes);
};

/// Split of the build cost, so the profile can attribute the intersection separately from the
/// canonicalisation. Filled only when a non-null pointer is passed.
struct BallGraphTiming {
    long long intersect_ns{0};
    long long finalize_ns{0};
};

/// Builds `H` for one shot.
///
/// `seeded_dets` must be the **post-preamble** detection events: the shot's events symmetric
/// -differenced with the graph's negative-weight detection events and with user-graph boundary
/// nodes removed, sorted ascending (§M2.1). `H` never sees a negative weight.
///
/// `horizon` is `T` in the flooder's time units — the same value handed to
/// `process_timeline_until_horizon`, per the unit rule of §0.
void build_ball_graph(
    const BallTables& tables,
    const std::vector<uint64_t>& seeded_dets,
    horizon_int horizon,
    BallGraphArena& arena,
    BallGraphBuildMode mode,
    BallGraphTiming* timing = nullptr);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_GRAPH_H
