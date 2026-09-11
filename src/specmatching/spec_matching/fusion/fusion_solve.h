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

/// §4 of the fusion design — the solver-side additions the divide-and-conquer solve needs, behind
/// one small interface.
///
/// The shape of the thing: a shot's whole `H` lives on **one** `pm::Mwpm`. An oversized component
/// has been cut into pieces of at most `n/k` defects by an edge cut; the cut edges are installed in
/// the instance but *masked*, and each cut defect carries a **dummy** boundary edge of weight
/// `floor(w_min/2)` standing in for them. Each piece is then injected into the shared instance at
/// whatever time the clock has reached and solved as an independent MWPM problem in which the cut is
/// a virtual boundary. A fusion unmasks the cut edges that have become internal, raises (or removes)
/// the dummy boundaries they were standing in for, releases the regions those dummies were holding
/// matched, and lets blossom continue from the combined primal-dual state.
///
/// Nothing here is a new algorithm. The two additions to the vendored solver are `GraphFlooder`'s
/// `edge_mask` and `dual_cap` (see `graph_flooder.h`); everything below is composition of the
/// flooder and matcher operations that already exist, which is what makes the result exact in
/// matching weight rather than approximately so.
///
/// **The horizon is a cap on duals, not on the clock.** In a monolithic run at horizon `T` the two
/// coincide, because every dual is at most the elapsed time. Here pieces are solved one after
/// another on one clock and released regions resume with radii they earned earlier, so the bound has
/// to live on each growing region. `attach` installs it; `TRUNCATED` means some defect's dual would
/// have passed `T`, which is exactly the condition under which edges of length `> 2T` — the ones `H`
/// does not contain — could become the next tight edge.

#ifndef SPECMATCHING_SPEC_MATCHING_FUSION_FUSION_SOLVE_H
#define SPECMATCHING_SPEC_MATCHING_FUSION_FUSION_SOLVE_H

#include <cstdint>
#include <limits>
#include <vector>

#include "specmatching/sparse_blossom/matcher/mwpm.h"
#include "specmatching/spec_matching/manifold/ball_graph.h"
#include "specmatching/spec_matching/truncation/horizon.h"
#include "specmatching/spec_matching/truncation/truncated_timeline.h"

namespace pm {
namespace spec_matching {
namespace fusion {

/// The dummy boundary weight standing in for a cut edge of weight `w`: the largest **even** integer
/// at most `floor(w / 2)`.
///
/// The design says `floor(w / 2)`, and the halving is what the soundness argument needs —
/// `bnd(u) + bnd(v) <= w` gives the crossing edge non-negative slack in the combined state, and
/// rounding down only ever helps that. The extra "even" is a parity rule the vendored solver
/// imposes and does not state:
///
///   Two regions growing towards each other along an edge of weight `w` are scheduled to meet at
///   `(w - y1 - y2) >> 1`. That shift is exact only when `w - d1 - d2` is even. Every edge weight in
///   the tree is even by construction (`horizon.h`: "every edge weight is even and every collision
///   happens at an integer time"), and in a monolithic run every dual starts at zero and parity
///   propagates along tight edges, so the condition holds everywhere it is used. It does **not**
///   hold across a fusion: two regions released from two pieces carry duals that two independent
///   solves chose, and if their sum is odd the shift rounds down and the pair is matched one tick
///   short of tight. The matching is then off by one unit of weight — which is exactly what a
///   `--verify` mismatch of 1 looks like.
///
///   A released region is always one that was matched to a dummy, and a boundary match is tight, so
///   its dual is exactly the dummy's weight. Making every dummy weight even therefore makes every
///   released dual even, and from an all-even start the solver's own propagation keeps it that way:
///   a tight edge `d_x + d_y = w` with `w` even forces `d_x` and `d_y` to agree in parity, so the
///   regions a fusion pulls in through stolen matches inherit it too.
///
/// `pm::weight_int` is unsigned, so the mask is the whole of the rounding.
inline pm::weight_int even_dummy_weight(pm::weight_int w) {
    return (w / 2) & ~(pm::weight_int)1;
}

/// The `boundary_weight` value that says "this defect has no boundary edge at all", i.e. the code
/// boundary is farther than `T` and no crossing edge is standing in for it.
constexpr pm::weight_int NO_BOUNDARY_WEIGHT = std::numeric_limits<pm::weight_int>::max();

/// §4.1's mask arena: one byte per *directed* adjacency entry of the shared instance, laid out CSR
/// style against `MatchingGraph::nodes`. Sized once at startup and re-offset per shot; nothing here
/// allocates once the high-water mark is reached.
struct EdgeMask {
    std::vector<uint8_t> bits;
    /// `offsets[i] .. offsets[i + 1]` is node `i`'s slice, parallel to `nodes[i].neighbors`.
    std::vector<uint32_t> offsets;
    /// Rewrites that had to enlarge a buffer, the mask-side twin of `BallGraphArena::grow_events`.
    uint64_t grow_events{0};

    void reserve(size_t n_max, size_t entries_max);

    /// Recomputes the offsets from the adjacency sizes the shot is *about* to write: `degrees[i]` is
    /// how many entries node `i` will end up with, boundary half-edge included.
    ///
    /// It has to be the degrees rather than the instance's own adjacency (which is what this used to
    /// read) because the graph is now built per piece inside the leaves — at the point the offsets
    /// are needed, before the first leaf runs, no node has any adjacency to measure. Nothing is
    /// cleared either: a per-piece build writes the mask byte of every entry it writes, so there is
    /// no bit left over from the previous shot to clear.
    void rewrite_offsets(const uint32_t* degrees, size_t num_nodes);

    inline uint32_t total_entries() const {
        return offsets.empty() ? 0 : offsets.back();
    }
    inline void set(uint32_t node, uint32_t neighbor_index, uint8_t value) {
        bits[offsets[node] + neighbor_index] = value;
    }
    inline uint8_t get(uint32_t node, uint32_t neighbor_index) const {
        return bits[offsets[node] + neighbor_index];
    }
    inline uint64_t capacity_signature() const {
        return (uint64_t)bits.capacity() + offsets.capacity() + grow_events;
    }
};

/// The shared instance, plus the state §4 needs to drive it. One of these per `(d, p)` cell; the
/// `pm::Mwpm` it points at is owned by the caller (a `BallMwpm`, rebuilt per shot).
struct FusionInstance {
    pm::Mwpm* mwpm{nullptr};
    EdgeMask mask;
    /// `T` in the flooder's time units.
    horizon_int horizon{0};

    /// The weight a boundary edge is raised to when its defect has no boundary left at all. A dual
    /// can never exceed `T`, so anything above `T` is unreachable, which is what "remove the
    /// boundary edge" means operationally; erasing neighbour 0 would be an `O(degree)` shift of four
    /// parallel arrays for no behavioural difference.
    ///
    /// `2T + 2` rather than `T + 1` for two reasons, both of which bite. It has to be **even**, for
    /// the parity rule in `even_dummy_weight` below. And it has to clear `T` by more than one tick:
    /// at `T + 1` the boundary collision would be scheduled for exactly the instant the dual cap
    /// fires, and which of the two the queue hands over first is not something this code should be
    /// deciding.
    inline pm::weight_int unreachable_boundary_weight() const {
        return (pm::weight_int)(2 * horizon + 2);
    }

    /// Points the flooder at the mask arena and installs the dual cap.
    ///
    /// Called **once at startup**, after the instance has been sized for `n_max` nodes, and not per
    /// shot. The per-shot attach this used to be was there because `BallMwpm::rebuild` replaces the
    /// whole `pm::Mwpm` when it has to grow the node pool, and with it the flooder holding these
    /// pointers — but with the graph built per piece (amendment 1 §3.1) there is no per-shot rebuild
    /// left to grow anything.
    void attach(pm::Mwpm& instance, size_t num_nodes, horizon_int T);

    /// Re-offsets the mask for a shot whose node `i` will have `degrees[i]` adjacency entries, and
    /// re-points the flooder at the arena.
    ///
    /// The re-pointing is the reason this is not just `mask.rewrite_offsets`: a shot with more
    /// adjacency entries than any before it reallocates `bits`, and the flooder is holding
    /// `bits.data()`. `attach` used to be what refreshed that pointer, and `attach` no longer runs
    /// per shot.
    void size_mask_for_shot(const uint32_t* degrees, size_t num_nodes);

#ifndef NDEBUG
    /// Amendment 1 test 7: one per adjacency record `build_piece` writes into a solver node — which
    /// is every solver-node write the *construction* of `H` makes, now that no whole-instance
    /// rebuild happens per shot. The profiler snapshots it around the manager's `scatter`, where it
    /// must not move.
    uint64_t debug_node_writes{0};
#endif
};

/// Amendment 1 §3.2 — everything one leaf's graph build reads, and the three things it writes.
///
/// Every pointer is into the caller's own fixed arenas; this owns nothing. The three written arrays
/// are partitioned by piece — a build touches only the entries belonging to its own defects — which
/// is what makes the builds order-independent (amendment 1 test 8).
struct PieceBuild {
    /// The shot's edge records, indexed by edge id, with `i < j`. Read in place: not copying these
    /// per shot is the point of the amendment.
    const BallGraphEdge* edges{nullptr};
    /// Each edge's observable mask. The design's input record is `(u, v, w, obs_mask)`;
    /// `BallGraphEdge` carries the ball-pool entry the mask is *derived* from, so the derivation
    /// happens once per shot with the rest of the input, outside every timed region, into this
    /// array. A build must not re-derive it from the ball tables.
    const pm::obs_int* edge_obs{nullptr};
    /// Piece id of every defect. An edge crosses a cut exactly when its endpoints disagree.
    const uint32_t* piece_of{nullptr};
    /// The piece-local boundary weight `min(w_b(u), floor(wmin(u) / 2))`, or `NO_BOUNDARY_WEIGHT`.
    const pm::weight_int* boundary_weight{nullptr};
    /// Whether that weight is a dummy standing in for a cut. A dummy carries no observable.
    const uint8_t* boundary_is_dummy{nullptr};
    /// The real code boundary's observable mask, used when the boundary is not a dummy.
    const pm::obs_int* boundary_obs{nullptr};
    /// `refused_slot[e]` is edge `e`'s index in the crossing-edge list, for the `cross_index_*`
    /// write below. Only read for an edge that crosses a cut.
    const uint32_t* refused_slot{nullptr};

    /// Written: the adjacency index of each crossing edge's two directions. The `lo` half belongs to
    /// the piece holding `edges[e].i` and the `hi` half to the piece holding `edges[e].j`, so two
    /// builds never write the same word.
    uint32_t* cross_index_lo{nullptr};
    uint32_t* cross_index_hi{nullptr};
    /// Written: one byte per defect, set once that defect's adjacency is complete.
    uint8_t* built{nullptr};
};

/// A solver core's scratch for `build_piece`. Its slices are indexed by the *mask's* offsets, so two
/// pieces can never touch the same word and one buffer per core is all a real system would need;
/// this binary runs one leaf at a time on one thread, so there is one of these.
struct PieceBuildScratch {
    /// The neighbour of each gathered adjacency entry, and the edge it came from.
    std::vector<uint32_t> neighbor;
    std::vector<uint32_t> edge;
    /// How many entries defect `u` has gathered so far. Only ever touched by `u`'s own piece.
    std::vector<uint32_t> cursor;

    void reserve(size_t n_max, size_t entries_max);
    /// Sizes the gather buffers to the shot's total adjacency entry count. Called from the manager's
    /// `scatter`, beside the mask it shares its offsets with.
    void size_for_shot(size_t entries, size_t num_nodes);

    inline uint64_t capacity_signature() const {
        return (uint64_t)neighbor.capacity() + edge.capacity() + cursor.capacity();
    }
};

/// Amendment 1 §3.2. Writes the adjacency, weights, observable masks, mask bytes and boundary edge
/// of every defect of piece `piece`, and of no other defect.
///
/// This is the work `BallMwpm::rebuild` used to do for the whole of `H` on the manager core, split
/// so that each leaf does its own part on its own core. Nothing here reads a value another piece's
/// build wrote: the pointer to a far node across a cut is an address into the shared node array and
/// is valid before that node is built, and the crossing edge is masked, so the flooder never follows
/// it before the fusion that unmasks it — by which time both sides are built.
///
/// `piece_edges` is the piece's slice of the edge list, in the non-decreasing weight order the input
/// arrives in. Adjacency has to come out in ascending neighbour index instead, because
/// `DetectorNode::index_of_neighbor` and the fusion's `neighbor_index` lookup are binary searches,
/// so each defect's gathered entries are ordered before they are written.
void build_piece(
    FusionInstance& fusion,
    PieceBuildScratch& scratch,
    const PieceBuild& in,
    uint32_t piece,
    const uint32_t* defects,
    size_t defect_count,
    const uint32_t* piece_edges,
    size_t piece_edge_count);

/// Moves the shared clock forward to the next even tick, if it is not already there.
///
/// The other half of `even_dummy_weight`'s parity rule, and it has the same one-line justification:
/// the flooder's growing/growing meeting point `(w - y1 - y2) >> 1` is exact only when the summed
/// duals are even. Every edge weight is even, so what has to hold is that a growing node's radius
/// agrees in parity with the clock — and that is preserved by every transition the solver makes
/// *given* that it holds when growth begins. Growth begins in exactly two places here: a piece is
/// injected (radius 0, so the clock must be even) and a fusion releases a region (radius equal to
/// its even dummy weight, so again the clock must be even). Between those points the queue is
/// empty — the previous solve returned only when it had drained — so the clock can be nudged
/// without disturbing a scheduled event.
///
/// Asserts the empty queue rather than checking for it: a caller that nudges the clock with events
/// outstanding would be scheduling into the past, which the monotone queue does not survive.
void align_clock_for_new_growth(pm::Mwpm& mwpm);

/// §4.3. Adds one detection event to the shared instance **at the current clock**, as a growing
/// region of zero radius with its dual cap installed.
///
/// `pm::begin_timeline` cannot be reused: it insists on an empty queue and resets the clock to zero,
/// and a piece injected after another piece has been solved has neither.
void inject_detection_event(pm::Mwpm& mwpm, horizon_int horizon, uint32_t node_index);

/// Runs the flooder until every region it can reach is matched, or until a dual cap fires.
///
/// Regions belonging to pieces that are already solved are all matched and frozen, and — with their
/// crossing edges masked — unreachable, so they are not touched.
TimelineStatus run_until_settled(pm::Mwpm& mwpm);

/// §4.3. `inject_detection_event` for each of `defects` (which must be ascending, so that the
/// equal-time tie-break inside the piece is the one it has in the full `H`), then
/// `run_until_settled`.
TimelineStatus solve_piece(pm::Mwpm& mwpm, horizon_int horizon, const uint32_t* defects, size_t count);

/// §4.4 step 1. Unmasks both directions of a crossing edge and reschedules its endpoints.
///
/// At the moment a fusion runs, every region of both children is matched or frozen, so a
/// frozen/frozen pair schedules nothing and this is bookkeeping. It is still done, because the
/// endpoints' next-event calculation is now allowed to see an edge it was not allowed to see
/// before, and the one place that must not be left to chance is an already-tight crossing edge:
/// once either endpoint is released below, its rescheduled scan finds the collision at the fused
/// timeline's current time. `bnd[u] + bnd[v] <= w` is what makes that "at or after now" rather than
/// "in the past".
void unmask_crossing_edge(
    FusionInstance& fusion, uint32_t u, uint32_t u_to_v_index, uint32_t v, uint32_t v_to_u_index);

/// §4.4 step 2. Sets defect `u`'s boundary edge to `(weight, obs_mask)` and, if `u`'s region is
/// matched to the boundary through `u`'s own boundary edge, **releases** it: the match is detached
/// and the region becomes a growing alternating-tree root with its radius unchanged.
///
/// Returns whether it released. A region holding two changed dummies releases once — the second call
/// finds it already in a tree and only rewrites the weight.
bool set_boundary_and_release(pm::Mwpm& mwpm, uint32_t u, pm::weight_int weight, pm::obs_int obs_mask);

/// §4.5. Extracts a finished tree node's matching: the observable mask and the total weight of the
/// pairs among `defects`.
///
/// Driven from the defect list rather than from the arena's live vector, because a shared instance
/// can hold several finished roots at once and `Harvester::extract_only_to_obs` would extract all of
/// them. Deduplication is free: extraction resets every node it owned, so a defect whose region has
/// already been extracted is reached with a null `region_that_arrived_top`.
pm::MatchingResult extract_root(pm::Mwpm& mwpm, const uint32_t* defects, size_t count);

/// §11 test 1, debug builds only: every unmasked edge of the instance has non-negative slack, and
/// every defect of `defects` has a dual of at most `T`.
///
/// Returns true when the invariant holds. Written as a predicate rather than an `assert` so the
/// caller can report *which* check failed and on which shot.
bool feasibility_holds(
    const FusionInstance& fusion, const uint32_t* defects, size_t count, const char** failure_out);

}  // namespace fusion
}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_FUSION_FUSION_SOLVE_H
