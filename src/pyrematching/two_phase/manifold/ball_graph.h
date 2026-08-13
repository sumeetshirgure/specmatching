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

/// The connected components of `H`, over its defect-defect edges only, and the undirected adjacency
/// the structural statistics walk.
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
    /// Members of component `c`, ascending: `[member_offsets[c], member_offsets[c + 1])` of
    /// `members`.
    std::vector<uint32_t> member_offsets;
    std::vector<uint32_t> members;
    /// Position of `H`-node `i` within its own component's block of `members`, so a per-component
    /// scratch array can be indexed `0..size-1` without a map.
    std::vector<uint32_t> local_index;

    /// Undirected CSR adjacency over `H`'s defect-defect edges, both directions. Boundary edges are
    /// not in it: the boundary is not a node, so it joins nothing (§A.2).
    std::vector<uint32_t> adj_offsets;
    std::vector<uint32_t> adj_target;
    std::vector<pm::weight_int> adj_weight;

    /// Per-component Dijkstra scratch for the weighted diameter, indexed by `local_index`.
    std::vector<pm::cumulative_time_int> dijkstra_dist;
    std::vector<uint8_t> dijkstra_done;
    /// Write cursors for the two counting-sort fills below (members, then adjacency). Reused rather
    /// than declared locally, so neither fill allocates in steady state.
    std::vector<uint32_t> fill_cursor;

    /// Analyses that had to grow a buffer, the component-side twin of `BallGraphArena::grow_events`.
    /// After a warmup shot this must stop increasing (invariant 7).
    uint64_t grow_events{0};

    inline size_t num_components() const {
        return roots.size();
    }
    inline uint32_t degree_of(uint32_t node) const {
        return adj_offsets[node + 1] - adj_offsets[node];
    }

    void clear();
};

/// Components larger than this are left without a weighted diameter and counted instead (§C.2). The
/// diameter is `O(s^2)` Dijkstras over the component, and at `p = 1e-3` a component this large is
/// already far outside the distribution the statistic is for.
inline constexpr uint32_t MAX_DIAMETER_COMPONENT_SIZE = 32;

/// §A's trivial-component prune: the union-find over `H`, the per-component verdict, and the
/// induced sub-`H` that the solver is actually handed.
///
/// This is **decode state**, not the profiling decomposition above. `BallComponents` builds the
/// full adjacency, the member blocks and the Dijkstra scratch that the §C statistics walk; this
/// builds only what §A.3 reads — a root per node, a size per root, and the one edge of every
/// two-member component — because it runs on the shot's own path. The two are deliberately not
/// shared: §C is computed for every shot whether or not the prune is enabled, and after the timed
/// window has closed, while this runs inside it.
///
/// No timer is started anywhere in the routines that fill this and none may be added: the prune is
/// assumed free (hardware-offloadable), and the latency account of this experiment is the solver
/// and the harvest running on a smaller node set.
struct BallPrune {
    /// Union-find parent while the unions run; the root of each `H`-node afterwards. Every link
    /// points at a **smaller** index, so a set's root is its minimum member and components
    /// enumerate in ascending root order without a sort (§0 determinism).
    std::vector<uint32_t> component_of;
    /// Size of the component rooted at an `H`-node. Only meaningful at a root, which is the only
    /// place it is read.
    std::vector<uint32_t> component_size;
    /// Index into `BallGraph::edges` of the single edge joining a two-member component, stored at
    /// that component's root; `NO_PAIR_EDGE` everywhere else. `H` holds each undirected pair once,
    /// so a two-member component has exactly one.
    std::vector<uint32_t> pair_edge;
    /// §A.3's verdict for the component rooted at an `H`-node, as a `TrivialVerdict`. Only
    /// meaningful at a root; the second pass spreads it to the members.
    std::vector<uint8_t> verdict;
    /// 1 for the `H`-nodes the solver still has to take (§A.4's SOLVER set).
    std::vector<uint8_t> node_to_solver;
    /// Index of an `H`-node in `solver_graph`, or `NOT_IN_SOLVER` when it was resolved away.
    std::vector<uint32_t> h_to_solver;
    /// The sub-graph of `H` induced on the SOLVER set, renumbered in ascending `H`-node order so
    /// that `h_to_det` stays strictly ascending and the residual maps back exactly as it does from
    /// the full `H` (§A.4).
    BallGraph solver_graph;

    /// Prunes that had to grow a buffer, the prune-side twin of `BallGraphArena::grow_events`.
    /// After a warmup shot this must stop increasing (invariant 7).
    uint64_t grow_events{0};

    static constexpr uint32_t NO_PAIR_EDGE = UINT32_MAX;
    static constexpr uint32_t NOT_IN_SOLVER = UINT32_MAX;

    void clear();
};

/// §A.2. Union-find over `H`'s defect-defect edges, and only those: the boundary is not a node of
/// `H`, so two defects that both reach it are not thereby connected, and joining them would merge
/// components the solver keeps apart.
///
/// Leaves `component_of` pointing every node straight at its root, `component_size` filled at the
/// roots, and `pair_edge` naming the one edge of each two-member component. `verdict` is sized and
/// zeroed for the caller to fill.
void compute_prune_components(const BallGraph& graph, BallPrune& prune);

/// §A.4. Fills `prune.solver_graph` with the sub-graph of `H` induced on the nodes flagged in
/// `prune.node_to_solver`, and `prune.h_to_solver` with the renumbering.
///
/// Both endpoints of an `H` edge are in one component and so share a verdict; the edge is kept iff
/// its component is. Boundary edges of kept nodes are carried over unchanged, which is what makes
/// the solver see exactly the problem it would have seen on the full `H` restricted to these
/// components.
void build_solver_subgraph(const BallGraph& graph, BallPrune& prune);

/// Owns everything `H`'s construction needs across shots, so that steady-state building allocates
/// nothing (debug invariant 11). `det_to_h` and the syndrome bitset are sized once to the detector
/// graph and cleared through touched lists, never by re-zeroing.
struct BallGraphArena {
    BallGraph graph;

    std::vector<uint32_t> det_to_h;
    std::vector<uint64_t> syndrome_words;
    std::vector<uint32_t> touched_words;

    /// §C's component decomposition of the shot's `H`. Filled only by `analyze_ball_components`,
    /// which the decode path never calls; it lives here for the arena's lifetime and reset
    /// discipline, so the analysis allocates nothing in steady state either (invariant 7).
    BallComponents components;

    /// §A's prune state and the sub-`H` it hands the solver. Unlike `components` this *is* on the
    /// decode path, and it lives here for the same reason: the arena's lifetime and reset
    /// discipline are what make a steady-state shot allocate nothing (invariant 11, extended).
    BallPrune prune;

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

/// Union-find over `H`'s defect-defect edges, then the members, the adjacency and the per-component
/// blocks the §C statistics read (§A.2).
///
/// Takes no timer and is given none: the component work is untimed by construction, and calling it
/// from inside a timed window would put it in a latency number it must never enter.
void analyze_ball_components(const BallGraph& graph, BallComponents& components);

/// The weighted diameter of component `index`: the largest shortest-path distance between two of its
/// members, with paths confined to the component's own subgraph of `H` and weighted by `w_int`.
///
/// This is an `H`-subgraph diameter, **not** a `G` diameter: a pair whose `G` geodesic leaves the
/// component is measured here by the route that stays inside it, and the two can differ. That is
/// fine for a structural metric — the quantity of interest is how far apart the component holds its
/// own members — but it is why the number must not be read as a distance in `G`.
///
/// Returns `-1` for a component larger than `MAX_DIAMETER_COMPONENT_SIZE`, which the caller counts
/// rather than computing (§C.2). A singleton has diameter 0.
///
/// `components` is taken by mutable reference for its Dijkstra scratch alone; nothing it describes
/// about `H` is modified.
pm::cumulative_time_int component_diameter(BallComponents& components, size_t index);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_GRAPH_H
