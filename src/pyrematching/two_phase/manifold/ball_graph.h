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

/// Element sizes of the ball-table pools the shot-time intersection reads, taken from the array
/// types themselves rather than assumed. The structural counters are quoted in bytes and the
/// conversion has to survive someone widening `pm::weight_int` or narrowing `ball_target`.
namespace ball_element_bytes {
inline constexpr uint64_t TARGET = sizeof(decltype(BallTables::ball_target)::value_type);
inline constexpr uint64_t WEIGHT = sizeof(decltype(BallTables::ball_w_int)::value_type);
inline constexpr uint64_t WORD = sizeof(decltype(BallTables::ball_words)::value_type);
inline constexpr uint64_t WORD_RANK = sizeof(decltype(BallTables::ball_word_rank)::value_type);
inline constexpr uint64_t ENTRY_BY_RANK = sizeof(decltype(BallTables::ball_entry_by_rank)::value_type);
inline constexpr uint64_t MASK_OFFSET = sizeof(decltype(BallTables::ball_mask_offsets)::value_type);
inline constexpr uint64_t MASK_ID = sizeof(decltype(BallTables::ball_mask_ids)::value_type);
inline constexpr uint64_t HAS_BCOST = sizeof(decltype(BallTables::has_bcost)::value_type);
inline constexpr uint64_t BCOST_WEIGHT = sizeof(decltype(BallTables::bcost_w_int)::value_type);
inline constexpr uint64_t BCOST_MASK_OFFSET = sizeof(decltype(BallTables::bcost_mask_offsets)::value_type);
inline constexpr uint64_t BCOST_MASK_ID = sizeof(decltype(BallTables::bcost_mask_ids)::value_type);
}  // namespace ball_element_bytes

/// What the intersection *moves*, rather than how long it takes on this laptop: the structural
/// (hardware-budget) counters of §M2. `intersect_ns` measures this machine's DRAM latency; these
/// measure the quantity a target architecture would have to move, and can be fed to a latency model
/// for any architecture. Filled only when a non-null pointer is passed.
///
/// The counters are accumulated at the sites that do the reads, so an early-terminating loop is
/// counted as it actually ran, not as the model would have it run.
struct BallGraphCounts {
    /// The traversal itself: the bitset window in `BITSET` mode, the ball entries walked in `SCAN`.
    /// Per-node CSR offset lookups (`ball_offsets`, `ball_word_offsets`, `ball_word_base`) are
    /// excluded in both modes — 16-20 B per defect against a ~1 kB window.
    uint64_t isect_scan_bytes{0};
    /// Fetched only because a candidate pair showed up: the rank -> entry indirection and the
    /// weight in `BITSET` mode (where the bitset carries neither), the observable id list of every
    /// pair that becomes an edge, and the `bcost_*` boundary lookups.
    uint64_t isect_hit_bytes{0};
    /// What the *other* build mode's traversal would have read on the same shot, derived from the
    /// tables without running it. This is the number that says whether the `SCAN`/`BITSET`
    /// crossover has moved at large `d`.
    uint64_t isect_scan_bytes_other_mode{0};
    /// Edge records emitted: undirected defect-defect pairs, each counted once, and boundary edges.
    uint64_t edges_written{0};
    uint64_t boundary_edges_written{0};

    inline uint64_t isect_bytes() const {
        return isect_scan_bytes + isect_hit_bytes;
    }
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
    BallGraphTiming* timing = nullptr,
    BallGraphCounts* counts = nullptr);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_GRAPH_H
