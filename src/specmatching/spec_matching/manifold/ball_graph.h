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

#ifndef SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_GRAPH_H
#define SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_GRAPH_H

#include <cstdint>
#include <vector>

#include "specmatching/spec_matching/manifold/ball_tables.h"

namespace pm {
namespace spec_matching {

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

/// The connected components of `H`, over its defect-defect edges only, and their sizes.
///
/// Size is the whole of it. The degrees, the two `H`-subgraph diameters and the boundary structure
/// this used to carry are gone, and with them the undirected adjacency, the per-component member
/// blocks and the Dijkstra/BFS scratch they were walked with — the only structural statistic the
/// `H` profiler reports is the component size distribution, so nothing else is built.
///
/// **Profiling state.** Nothing on the decode path reads any of this, no timer is ever started
/// inside the routine that fills it, and it is never filled from inside a timed window: the
/// component work is assumed free (hardware-offloadable), so charging a latency number for it would
/// be measuring a stage this experiment does not intend to run on a CPU.
///
/// The root of a component is the **smallest** `H`-node index in it, so components enumerate in
/// ascending root order and every list below is a function of `H` alone (§0 determinism).
struct BallComponents {
    /// Root of `i`'s component, for each `H`-node `i`. Doubles as the union-find parent array while
    /// the unions are running.
    std::vector<uint32_t> component_of;
    /// Index of `i`'s component within `roots` / `sizes` / `member_offsets`, for each `H`-node `i`.
    /// Only meaningful at a root, which is where the fill reads it.
    std::vector<uint32_t> component_index;
    /// Roots, ascending. One entry per component.
    std::vector<uint32_t> roots;
    std::vector<uint32_t> sizes;

    /// Analyses that had to grow a buffer, the component-side twin of `BallGraphArena::grow_events`.
    /// After a warmup shot this must stop increasing (invariant 7).
    uint64_t grow_events{0};

    inline size_t num_components() const {
        return roots.size();
    }

    void clear();
};

/// §1/§2. The connected components of `H`, and the per-component sub-`H` each one is solved on.
///
/// This is **decode state**: the solver no longer sees `H`, it sees one component at a time. What
/// licenses that is the independence property of §1 — two regions interact only when the sum of
/// their radii reaches the distance between their defects, every radius is `<= T` by the horizon,
/// so interaction needs `d_G <= 2T`, which is exactly the condition for an `H` edge. Defects in
/// different components have no `H` edge and so never interact before `T`; the boundary is not a
/// node, does not grow, and a boundary match in one component changes no state in any other.
///
/// It is deliberately **not** `BallComponents`, which is the statistics' own decomposition and
/// carries nothing but sizes. This builds what the solve needs — member blocks, internal-edge blocks
/// and the boundary lookup — and the statistics join to it by component index, which is well-defined
/// because both enumerate in ascending root order.
struct BallComponentSplit {
    /// Union-find parent while the unions run; the root of each `H`-node afterwards. Every link
    /// points at a **smaller** index, so a set's root is its minimum member and components
    /// enumerate in ascending root order without a sort (§0 determinism).
    std::vector<uint32_t> component_of;
    /// Index of `i`'s component within `roots` / `sizes` / `member_offsets`. Only meaningful at a
    /// root, which is where the fills read it.
    std::vector<uint32_t> component_index;
    /// Roots, ascending. One entry per component.
    std::vector<uint32_t> roots;
    std::vector<uint32_t> sizes;
    /// Members of component `c`, ascending: `[member_offsets[c], member_offsets[c + 1])`. Ascending
    /// is the whole of §1's bit-exactness clause — relabelling each component's nodes in the same
    /// ascending order they have in the full `H` leaves blossom's equal-time tie-break unchanged.
    std::vector<uint32_t> member_offsets;
    std::vector<uint32_t> members;
    /// Position of `H`-node `i` within its own component's member block, i.e. its index in the
    /// sub-`H`.
    std::vector<uint32_t> local_index;

    /// Indices into `BallGraph::edges` of the edges internal to component `c`, ascending:
    /// `[edge_offsets[c], edge_offsets[c + 1])`. Both endpoints of an `H` edge are in one
    /// component, so every edge appears in exactly one block.
    std::vector<uint32_t> edge_offsets;
    std::vector<uint32_t> edge_index;

    /// Index into `BallGraph::boundary_edges` of `H`-node `i`'s boundary edge, or `NO_BOUNDARY`.
    /// `H` gives a node at most one, so this is a lookup rather than a range: §2.4's rule is to
    /// copy `H`'s boundary edges for the component's members and not to re-derive the cutoff.
    std::vector<uint32_t> node_boundary;

    /// `TimelineStatus` of component `c`, as `COMPONENT_COMPLETE` / `COMPONENT_TRUNCATED`. Written
    /// by the decode, read by the escalation predicate and by the §3.5 statistics. Kept as a byte
    /// so that this header does not have to reach into the matcher for an enum.
    std::vector<uint8_t> status;

    /// The sub-`H` of the component being solved right now, rewritten in place per component.
    BallGraph sub;

    /// Write cursors for the two counting-sort fills (members, then internal edges). Reused rather
    /// than declared locally, so neither fill allocates in steady state.
    std::vector<uint32_t> fill_cursor;

    /// Decompositions that had to grow a buffer, the split-side twin of
    /// `BallGraphArena::grow_events`. After a warmup shot this must stop increasing (invariant 11).
    uint64_t grow_events{0};

    static constexpr uint32_t NO_BOUNDARY = UINT32_MAX;
    static constexpr uint8_t COMPONENT_COMPLETE = 0;
    static constexpr uint8_t COMPONENT_TRUNCATED = 1;

    inline size_t num_components() const {
        return roots.size();
    }
    inline bool any_truncated() const {
        for (uint8_t s : status) {
            if (s == COMPONENT_TRUNCATED)
                return true;
        }
        return false;
    }
    inline size_t num_truncated() const {
        size_t count = 0;
        for (uint8_t s : status)
            count += s == COMPONENT_TRUNCATED ? 1 : 0;
        return count;
    }

    void clear();
};

/// §2.2. Union-find over `H`'s defect-defect edges, and only those: the boundary is not a node of
/// `H`, so two defects that both reach it are not thereby connected, and joining them would merge
/// components the solve keeps apart.
///
/// Leaves the components enumerated in ascending root order, the member and internal-edge blocks
/// filled ascending, the boundary lookup filled, and `status` sized and zeroed for the decode to
/// write.
void decompose_ball_components(const BallGraph& graph, BallComponentSplit& split);

/// §2.4. Rewrites `split.sub` as the sub-graph of `H` induced on component `index`.
///
/// Local node `j` is `members[member_offsets[index] + j]`, so `sub.h_to_det` is a strictly
/// ascending subsequence of `graph.h_to_det` — which is what keeps the tie-break order within the
/// component identical to the one it has in the full `H`, and what lets a residual sorted in the
/// sub-graph be sorted in `G`. Edge weights, observable entries and boundary cutoffs are **copied**
/// from `H`; nothing is re-derived here.
void build_component_subgraph(const BallGraph& graph, BallComponentSplit& split, size_t index);

/// Owns everything `H`'s construction needs across shots, so that steady-state building allocates
/// nothing (debug invariant 11). `det_to_h` and the syndrome bitset are sized once to the detector
/// graph and cleared through touched lists, never by re-zeroing.
struct BallGraphArena {
    BallGraph graph;

    std::vector<uint32_t> det_to_h;
    std::vector<uint64_t> syndrome_words;
    std::vector<uint32_t> touched_words;

    /// §3.5.2's component decomposition of the shot's `H` — the adjacency, the diameters and the
    /// distributions. Filled only by `analyze_ball_components`, which the decode path never calls;
    /// it lives here for the arena's lifetime and reset discipline, so the analysis allocates
    /// nothing in steady state either (invariant 7).
    BallComponents components;

    /// §2's component split and the sub-`H` the solve runs on. Unlike `components` this *is* the
    /// decode path, and it lives here for the same reason: the arena's lifetime and reset
    /// discipline are what make a steady-state shot allocate nothing (invariant 11, extended).
    BallComponentSplit split;

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

/// Union-find over `H`'s defect-defect edges, then the components in ascending root order and their
/// sizes. That is the whole of §3.5.2 now: the size of each component, and nothing else about it.
///
/// Takes no timer and is given none: the component work is untimed by construction, and calling it
/// from inside a timed window would put it in a latency number it must never enter.
void analyze_ball_components(const BallGraph& graph, BallComponents& components);

}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_MANIFOLD_BALL_GRAPH_H
