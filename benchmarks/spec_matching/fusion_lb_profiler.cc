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

/// `fusion_lb_profiler` — divide-and-conquer sparse blossom on `H` with a `k`-core list scheduler.
///
/// The load balancer priced a shot whose components were spread over `k` cores; its ceiling is the
/// largest component, which scheduling cannot split. This binary splits it. Following the
/// division/fusion idea of Fusion Blossom (Wu & Zhong, 2023), an oversized component is cut into
/// pieces of at most `n/k` defects, each piece is solved as an independent MWPM problem in which the
/// cut is a virtual boundary, and the pieces are fused back together by continuing blossom from
/// their combined primal-dual state. The cut is an edge cut of `H`, chosen inside the union-find
/// pass that already finds the components, and it exists only in software: no ball-table path
/// changes. The result is exact in matching weight.
///
/// ```
///   sparse_k = makespan of the simulated schedule (preprocessing, leaves, fusions, extraction, combine)
///   system   = escalated ? fallback : min(fallback, sparse_k)
/// ```
///
/// `fallback` is the stock decoder on the full detector graph `G`, timed on its own core outside
/// the `k`.
///
/// ## Execution and memory model — `fusion_k_core_model`
///
/// **One thread; no thread is ever spawned.** Every unit of work runs on the profiling thread and is
/// timed with `hires_timer.h`'s thread-CPU clock; a simulated list scheduler places those
/// measurements on `k` solver cores plus one resource-manager core plus one fallback core, and the
/// per-shot figure is the makespan of that schedule.
///
/// **One shared solver instance per shot.** A shot's whole `H` — every region, every event queue —
/// lives in one `BallMwpm` that every simulated core addresses directly. Nothing is copied between
/// cores and the transfer cost is modelled as **zero**; the cache-miss cost a real core would pay to
/// touch another core's data is not modelled. The warm-cache caveat travels with every row: a unit
/// of work executed immediately after another sees a warmer cache than a real core would.
///
/// **Nothing outside a timed region is on the critical path** — the ball-table edge list, the
/// instance reset, the output writing and the `--verify` cross-check are all outside.
///
/// ## Deviations from the design, recorded rather than silently taken
///
///  1. **Location and corpus.** The file is `benchmarks/spec_matching/fusion_lb_profiler.cc`, not
///     the design's `benchmarks/two_phase/`, which is this tree's name for that corpus; the
///     generator call is the corpus's own (`profiler_util.h`'s `ShotSampler::make`), and it carries
///     that call's two inherited deviations — the basis is `rotated_memory_x`, and three noise
///     parameters are set rather than four. `--T` is in **multiples of one lattice edge weight**,
///     as in every other binary in this directory, not in raw DEM float units; `run.log` records
///     `T`, `T_weight_units` and `T_int`.
///  2. **The edge list arrives sorted (amendment 1 §1).** The hardware that produces the edge list
///     sorts it, so the solver receives `H`'s edges in non-decreasing weight, ties by `(u, v)`, and
///     no bucketing pass exists. The simulator stands in for that sort with a `std::sort` on the
///     edge indices, run **outside every timed region** with the rest of the input; `ticks_bucket`
///     is gone from every output. The refusal property of §3.2 needs only non-decreasing order,
///     which this provides, and `(w, u, v)` is a strict total order on `H`'s edges, so the
///     permutation is the one the design's bucketing would have produced.
///  3. **When the quotient graph is built (§3.2/§3.3).** The design records a refused edge against
///     the pair of union-find roots its endpoints have *at the moment of refusal*. Those roots are
///     not final — a refusal between `a` and `b` does not stop `a` from later merging with some
///     smaller `c` — so keying on them would key on a stale name. The refused edges are recorded as
///     edge indices during the pass (which is all the `bnd` update needs, and that update is
///     unchanged) and the quotient graph is built from their **final** piece ids immediately
///     afterwards, inside the same timed `uf` region. Nothing is lost: two pieces separated by a
///     refusal can never merge later, because set sizes only grow and the refusal means their sum
///     already exceeded `S`.
///  4. **When a fusion's interior edge list is built (§3.3).** The design has each internal node
///     accumulate its list as the greedy merges adjacency. Keeping one list per quotient edge
///     correct under merging needs it stored on exactly one side of each adjacency pair, which the
///     small-to-large merge does not give you. The greedy therefore merges **counts** only, and each
///     refused edge is assigned to the internal node that is the lowest common ancestor of its two
///     pieces' leaves once the tree is finished. That node is by definition the fusion at which the
///     edge becomes interior, so the lists are the same lists.
///  5. **An escalating shot stops where it escalated (§5).** The design runs the remaining jobs
///     after a `TRUNCATED` one "so per-node numbers are complete". A truncated node leaves
///     alternating trees standing and a primal-dual state that §4.4 is not defined on, so fusing
///     into it is not a measurement of anything a decoder would do. The shot is abandoned at the
///     first truncation instead, the row records how far it got, and `agg.json` takes its `ticks_*`
///     and `cp_*` means over **non-escalating shots only**. `sparse_k` was never going to be used on
///     an escalating shot in any case.
///  6. **Per-node adjacency ordering (amendment 1 §3.2).** A piece's edge list arrives in weight
///     order and adjacency has to be written in ascending neighbour index, so each defect's entries
///     are gathered and then **insertion-sorted by neighbour** inside the leaf — the design's own
///     "acceptable and simpler" option, on degrees that are small. `run.log` records the choice.
///     The gather slices are the edge mask's own offsets, so no two pieces address the same word
///     and the scratch needs no per-piece bookkeeping at all.
///  7. **"Removing" a boundary edge (§4.4 step 2).** A defect that ends a fusion with neither a real
///     boundary nor a crossing edge left has its boundary half-edge raised to `2T + 2` rather than
///     erased. A dual can never exceed `T`, so the edge is unreachable, which is what removal means
///     operationally; erasing neighbour 0 would shift four parallel arrays for no behavioural
///     difference.
///  8. **The dummy boundary weight is `floor(w / 2)` rounded down to even.** The design says
///     `floor(w / 2)`; the extra rounding is a parity rule the vendored solver imposes and does not
///     state, and without it `--verify` fails by exactly one unit of weight on a small fraction of
///     shots. Two regions growing towards each other meet at `(w - y1 - y2) >> 1`, which is exact
///     only when the summed duals are even; a monolithic run gets that for free (all duals start at
///     zero, all edge weights are even, and parity propagates along tight edges), while two pieces
///     solved independently do not. A released region was matched to a dummy and a boundary match is
///     tight, so its dual *is* the dummy's weight — making every dummy even makes every released
///     dual even, and the solver's own propagation carries it from there. Rounding down only
///     strengthens the soundness inequality `bnd(u) + bnd(v) <= w`. See
///     `fusion::even_dummy_weight`, and `fusion::align_clock_for_new_growth` for the other half of
///     the same rule: growth also has to *begin* on an even clock, which it does because the queue is
///     drained between one solve and the next and the clock can be nudged there for free.
///  9. **§8's escalation clause is counted, not asserted.** The design asks `--verify` to check that
///     the fused solve and the monolithic one "both escalate or neither does". That is not a
///     theorem, and it fails on real shots: truncation at `T` is a property of the *trajectory* the
///     solver takes, not of the problem, and a fused solve takes a different trajectory — its pieces
///     reach a dual optimum of their own, and a released region resumes from a radius the monolithic
///     run never gave it. So `--verify` gates on the **weight identity**, which is the exactness
///     claim, and reports the escalation difference as four counters in `agg.json`
///     (`escalation_divergence_*`).
///
///     Both directions are safe, which is why this is a rate and not a bug. A shot the fused solve
///     escalates goes to the fallback on `G`, which is exact. A shot it completes carries its own
///     certificate — a feasible dual for `H` with every dual at most `T`, and a perfect matching on
///     tight edges — and "every dual at most `T`" is exactly what makes `H`'s `2T` edge set
///     sufficient, so that answer is a global minimum-weight matching however the monolithic
///     trajectory went. What the counters measure is the escalation *rate*, which is a latency
///     question and belongs in a campaign, not in a gate.
///
/// ## What is timed
///
/// ```
///     compute_seeded_detection_events + build_ball_graph        the input       NOT timed
///     the hardware sort of the edge list, and the per-edge observable masks       NOT timed
///     marking every node of the shared instance unbuilt                           NOT timed
///   [ PRE ] on the manager core, three sub-regions sharing their boundary reads:
///     [ uf ]       bounded union-find, boundary weights, quotient graph
///     [ tree ]     cheapest-boundary-first fusion tree, then the interior edge lists
///     [ scatter ]  the sorted edge list grouped by piece, and the bookkeeping a fusion reads.
///                  No ball table, no solver node: O(n + m) integer work in fixed arenas.
///   [ DISPATCH ]  on the manager core, once per assigned job: the ready-heap pop and the core pick
///   [ LEAF(i) ]   build_piece(i) then solve_piece(i), two regions   on one solver core
///   [ FUSE(p) ]   fuse(p)                                                   on a solver core
///   [ EXTRACT(r) ] extract(r)                              on r's core, no dispatch charged
///   [ COMBINE ]   XOR of root observable masks, sum of root weights         on the manager core
///   [ FALLBACK ]  pm::decode_detection_events on G from the full syndrome   on the fallback core
///     the monolithic cross-check under --verify                                    NOT timed
///     abandon_shot / reset_for_next_shot                                           NOT timed
/// ```
///
/// The graph build is **in the leaf**, not on the manager (amendment 1). The manager used to rebuild
/// the whole of `H` — every node's adjacency, weights and observable masks for the entire shot,
/// serially — inside `scatter`, which put a stage that does not shrink with `k` on the critical
/// path. Each piece now writes its own defects' adjacency on its own core, and the shared instance
/// is sized and attached once at startup rather than rebuilt per shot.
///
/// A warm-up shot runs the identical workload and reads no clock at all: `run_shot` is templated on
/// the switch, so the reads compile out of the warm-up instantiation. Rows are stored
/// **uncorrected**; `timer_overhead_ns` goes to `run.log` and the analysis subtracts it.
///
/// ## Usage
///
///   fusion_lb_profiler [--d 5,7,9] [--p 1e-3] [--T 1.5,2] [--k 2,4,8] [--alpha 1]
///                      [--shots N] [--warmup 1000] [--seed S] [--verify]
///                      [--check-determinism] [--out DIR]
///
/// Writes one `shots_d{d}_p{p}_T{T}_k{k}.csv` per cell, plus `agg.json` and `run.log`, to `--out`.
/// The agent runs no campaigns; smoke runs at `--shots 100 --d 5` exist only to check that the
/// output parses, and their numbers appear in no document.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmarks/spec_matching/hires_timer.h"
#include "benchmarks/spec_matching/profiler_util.h"
#include "specmatching/spec_matching/fusion/fusion_solve.h"
#include "specmatching/spec_matching/manifold/ball_decoding.h"

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#ifndef SPECMATCHING_GIT_HASH
#define SPECMATCHING_GIT_HASH "unknown"
#endif
#ifndef SPECMATCHING_BUILD_FLAGS
#define SPECMATCHING_BUILD_FLAGS "unknown"
#endif
#ifndef SPECMATCHING_STIM_VERSION
#define SPECMATCHING_STIM_VERSION "unknown"
#endif

using namespace pm::spec_matching;
using namespace pm::spec_matching::profiler;
namespace fus = pm::spec_matching::fusion;

namespace {

/// Shots sampled per batch, as in the rest of the corpus. The stream is a function of the seed and
/// the batch sizes alone, and both are identical across the `T` and `k` sweeps, which is what makes
/// shot `i` of `(d, p)` the same shot in every cell.
constexpr size_t SHOT_BATCH = 4096;

constexpr size_t ROW_ARENA_CAPACITY = 1u << 17;
constexpr size_t ROW_ARENA_FLUSH_AT = ROW_ARENA_CAPACITY / 2;

/// Starting capacity for everything sized by `|E(H)|`, which has no static bound. The warm-up grows
/// these to the corpus's high-water mark exactly as `BallGraphArena` and `BallMwpm` do, and
/// `run.log` reports `no_growth_after_warmup` per cell.
constexpr size_t INITIAL_EDGE_CAPACITY = 4096;

/// "No boundary edge at all": the code boundary is farther than `T` and no crossing edge is standing
/// in for it. The same value the build reads it as.
constexpr pm::weight_int WEIGHT_INF = fus::NO_BOUNDARY_WEIGHT;
constexpr uint32_t NONE = UINT32_MAX;

struct Options {
    std::vector<size_t> distances = {5};
    std::vector<double> error_rates = {1e-3};
    std::vector<double> horizons = {1.5, 2.0};
    std::vector<size_t> core_counts = {2, 4, 8};
    /// Leaf cost exponent for the scheduler's priority. At `alpha == 1` the leaf cost is the defect
    /// count and no floating-point operation happens anywhere on the critical path; otherwise a
    /// table `cost[s] = llround(s^alpha)` is filled once at startup and looked up.
    double alpha = 1.0;
    size_t shots = 10000;
    size_t warmup = 1000;
    uint64_t seed = 20260910;
    bool verify = false;
    bool check_determinism = false;
    std::string out_dir = "benchmarks/spec_matching/results/fusion_lb_profiler";
};

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--d" || flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--p" || flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--T" || flag == "--horizons") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--k" || flag == "--cores") {
            options.core_counts = parse_list<size_t>(next());
        } else if (flag == "--alpha") {
            options.alpha = std::stod(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--warmup") {
            options.warmup = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--verify") {
            options.verify = true;
        } else if (flag == "--check-determinism") {
            options.check_determinism = true;
        } else if (flag == "--out" || flag == "--out-dir") {
            options.out_dir = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    if (options.distances.empty() || options.error_rates.empty() || options.horizons.empty() ||
        options.core_counts.empty())
        throw std::invalid_argument("--d, --p, --T and --k each need at least one value");
    if (options.shots == 0)
        throw std::invalid_argument("--shots must be positive");
    for (size_t distance : options.distances) {
        if (distance < 3 || distance % 2 == 0)
            throw std::invalid_argument("--d values must be odd and at least 3");
    }
    for (double noise : options.error_rates) {
        if (!(noise > 0) || !(noise < 1))
            throw std::invalid_argument("--p values must lie in (0, 1)");
    }
    for (double horizon : options.horizons) {
        if (!(horizon > 0))
            throw std::invalid_argument("--T values must be positive");
    }
    for (size_t k : options.core_counts) {
        if (k == 0)
            throw std::invalid_argument("--k values must be positive");
    }
    if (!(options.alpha > 0))
        throw std::invalid_argument("--alpha must be positive");
    std::sort(options.core_counts.begin(), options.core_counts.end());
    options.core_counts.erase(
        std::unique(options.core_counts.begin(), options.core_counts.end()), options.core_counts.end());
    return options;
}

std::string fmt_g(double value) {
    char buffer[40];
    std::snprintf(buffer, sizeof(buffer), "%.10g", value);
    return std::string(buffer);
}

/// `__int128` has no stream inserter, and a sum of squares of nanosecond ticks over a campaign
/// leaves the range a JSON number is guaranteed to survive, so the 128-bit accumulators print as
/// JSON *strings*.
std::string fmt_i128(__int128 value) {
    if (value == 0)
        return "0";
    bool negative = value < 0;
    unsigned __int128 magnitude = negative ? (unsigned __int128)(-value) : (unsigned __int128)value;
    char buffer[48];
    size_t at = sizeof(buffer);
    while (magnitude != 0) {
        buffer[--at] = (char)('0' + (int)(magnitude % 10));
        magnitude /= 10;
    }
    std::string out = negative ? "-" : "";
    out.append(buffer + at, sizeof(buffer) - at);
    return out;
}

/// FNV-1a over the deterministic columns of a cell (§10). Two runs of the same
/// `(seed, d, p, T, k, alpha)` must agree on it.
struct Digest {
    uint64_t value{1469598103934665603ull};

    void bytes(const void* pointer, size_t count) {
        const uint8_t* at = (const uint8_t*)pointer;
        for (size_t i = 0; i < count; i++) {
            value ^= at[i];
            value *= 1099511628211ull;
        }
    }
    template <typename T>
    void add(const T& item) {
        bytes(&item, sizeof(T));
    }
};

// ---------------------------------------------------------------------------------------------
// §3. The manager's preprocessing.
// ---------------------------------------------------------------------------------------------

/// One entry of a tree node's quotient adjacency during the §3.3 greedy. `other` is a tree-node id
/// that may since have been merged away; it is resolved through `tuf` at use, which is what lets the
/// merge touch only the node it is building and not its neighbours' lists.
struct QuotientAdj {
    uint32_t other;
    uint32_t count;
};

/// A candidate merge on the greedy's min-heap, keyed on `(count, min_id, max_id)`.
struct MergeCandidate {
    uint32_t count;
    uint32_t lo;
    uint32_t hi;

    /// `std::push_heap` builds a max-heap, so this is reversed: `a` orders *after* `b` when `a` is
    /// the cheaper merge, which puts the cheapest at the top.
    bool operator<(const MergeCandidate& other) const {
        if (count != other.count)
            return count > other.count;
        if (lo != other.lo)
            return lo > other.lo;
        return hi > other.hi;
    }
};

/// §3's fixed arenas: the sorted edge order, the bounded union-find, the quotient graph, the fusion
/// tree and the per-piece edge lists, all sized at cell setup and reused for every shot.
struct Cut {
    // ---- amendment 1 §1: the edge list in non-decreasing weight, as the hardware hands it over.
    // Filled outside every timed region; no bucket arena exists.
    std::vector<uint32_t> order;

    // ---- §3.2
    std::vector<uint32_t> parent;
    std::vector<uint32_t> set_size;
    /// `H`'s own boundary edge for each defect, or `WEIGHT_INF` when the code boundary is farther
    /// than `T`, plus its observable mask.
    std::vector<pm::weight_int> w_boundary;
    std::vector<pm::obs_int> obs_boundary;
    /// The piece-local boundary weight, `min(w_b(u), floor(wmin(u) / 2))`, and whether the second
    /// term won (i.e. whether the boundary edge currently in the instance is a dummy).
    std::vector<pm::weight_int> bnd;
    std::vector<uint8_t> via_dummy;
    std::vector<uint32_t> refused;

    // ---- pieces
    std::vector<uint32_t> piece_of;
    std::vector<uint32_t> piece_index_of_root;
    std::vector<uint64_t> root_stamp;
    uint64_t epoch{0};
    std::vector<uint32_t> piece_offsets;
    std::vector<uint32_t> piece_members;
    std::vector<uint32_t> fill_cursor;
    uint32_t n_pieces{0};
    uint32_t largest_piece{0};

    // ---- crossing edges incident to each defect (CSR)
    std::vector<uint32_t> cross_offsets;
    std::vector<uint32_t> cross_edges;

    // ---- amendment 1 §2.2: the sorted edge list grouped by piece, so that a leaf can build its own
    // part of the graph without scanning the whole list. An interior edge appears once, under its
    // endpoints' one piece; a crossing edge appears twice, once under each endpoint's piece, because
    // each endpoint writes its own adjacency entry and only its own.
    std::vector<uint32_t> piece_edge_offsets;
    std::vector<uint32_t> piece_edges;
    /// How many adjacency entries each defect will end up with, boundary half-edge included. The
    /// edge mask and the build scratch are both laid out against this.
    std::vector<uint32_t> node_degree;

    // ---- §3.3 quotient graph, open addressing on the piece pair
    std::vector<uint64_t> qkey;
    std::vector<uint32_t> qcount;
    std::vector<uint32_t> qslots;
    size_t qcap{0};

    // ---- §3.3 fusion tree
    std::vector<uint32_t> tree_child0;
    std::vector<uint32_t> tree_child1;
    std::vector<uint32_t> tree_parent;
    std::vector<uint32_t> tree_cost;
    std::vector<uint32_t> tree_depth;
    std::vector<uint8_t> tree_alive;
    std::vector<uint32_t> tree_roots;
    uint32_t n_tree_nodes{0};
    uint32_t n_fusions{0};
    uint32_t max_fusion_boundary{0};
    uint32_t depth{0};
    /// Union-find over tree nodes. Used twice: by the greedy, to resolve a stale adjacency entry to
    /// the live node it has been merged into, and again at run time, to answer "is this defect
    /// inside the node being fused?". Reset between the two.
    std::vector<uint32_t> tuf;
    std::vector<std::vector<QuotientAdj>> adj;
    std::vector<uint32_t> adj_touched;
    std::vector<uint64_t> mark_stamp;
    std::vector<uint32_t> mark_pos;
    uint64_t mark_epoch{0};
    std::vector<MergeCandidate> heap;

    // ---- the crossing edges each fusion turns interior (CSR over internal nodes)
    std::vector<uint32_t> interior_offsets;
    std::vector<uint32_t> interior_edges;

    // ---- §5 scheduling inputs
    std::vector<uint64_t> cost_tab;

    /// Scratch for the LCA walk of deviation 4, and for the per-defect dedup inside a fusion.
    std::vector<uint64_t> defect_stamp;
    uint64_t defect_epoch{0};

    void configure(size_t n_max, double alpha) {
        parent.assign(n_max, 0);
        set_size.assign(n_max, 0);
        w_boundary.assign(n_max, WEIGHT_INF);
        obs_boundary.assign(n_max, 0);
        bnd.assign(n_max, WEIGHT_INF);
        via_dummy.assign(n_max, 0);
        piece_of.assign(n_max, 0);
        piece_index_of_root.assign(n_max, 0);
        root_stamp.assign(n_max, 0);
        piece_offsets.assign(n_max + 1, 0);
        piece_members.assign(n_max, 0);
        // Indexed by *tree* node in `build_tree` as well as by piece and by defect, so it is sized
        // for the largest of the three.
        fill_cursor.assign(2 * n_max + 2, 0);
        cross_offsets.assign(n_max + 1, 0);
        piece_edge_offsets.assign(n_max + 1, 0);
        node_degree.assign(n_max, 0);
        defect_stamp.assign(n_max, 0);
        epoch = 0;
        defect_epoch = 0;

        size_t tree_max = 2 * n_max + 1;
        tree_child0.assign(tree_max, NONE);
        tree_child1.assign(tree_max, NONE);
        tree_parent.assign(tree_max, NONE);
        tree_cost.assign(tree_max, 0);
        tree_depth.assign(tree_max, 0);
        tree_alive.assign(tree_max, 0);
        tuf.assign(tree_max, 0);
        adj.assign(tree_max, {});
        mark_stamp.assign(tree_max, 0);
        mark_pos.assign(tree_max, 0);
        interior_offsets.assign(tree_max + 1, 0);
        tree_roots.reserve(n_max);
        adj_touched.reserve(tree_max);
        heap.reserve(INITIAL_EDGE_CAPACITY);

        order.reserve(INITIAL_EDGE_CAPACITY);
        refused.reserve(INITIAL_EDGE_CAPACITY);
        qslots.reserve(INITIAL_EDGE_CAPACITY);
        cross_edges.reserve(2 * INITIAL_EDGE_CAPACITY);
        piece_edges.reserve(2 * INITIAL_EDGE_CAPACITY);
        interior_edges.reserve(INITIAL_EDGE_CAPACITY);

        // Sized here rather than left to the first shot that needs it, so that the quotient table is
        // not the thing a cell's `no_growth_after_warmup` trips over: it doubles on the first shot
        // with more refused edges than any before, which for a heavy-tailed corpus can be well past
        // the warm-up.
        size_quotient_table(INITIAL_EDGE_CAPACITY);

        cost_tab.assign(n_max + 1, 0);
        for (size_t s = 0; s <= n_max; s++)
            cost_tab[s] = alpha == 1.0 ? (uint64_t)s : (uint64_t)std::llround(std::pow((double)s, alpha));
    }

    /// Sized to the shot's edge count; the hash never rehashes, and a full table is impossible
    /// because the number of distinct piece pairs is bounded by the number of refused edges.
    void size_quotient_table(size_t m) {
        size_t want = 4;
        while (want < 4 * (m + 1))
            want <<= 1;
        if (want > qcap) {
            qcap = want;
            qkey.assign(qcap, UINT64_MAX);
            qcount.assign(qcap, 0);
        }
    }

    inline uint32_t find(uint32_t node) {
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        return node;
    }
    inline uint32_t find_tree(uint32_t node) {
        while (tuf[node] != node) {
            tuf[node] = tuf[tuf[node]];
            node = tuf[node];
        }
        return node;
    }

    uint64_t capacity_signature() const {
        return (uint64_t)order.capacity() + refused.capacity() + cross_edges.capacity() +
               piece_edges.capacity() + interior_edges.capacity() + heap.capacity() + qkey.capacity() +
               tree_roots.capacity() + adj_touched.capacity() + qslots.capacity();
    }
};

// ---------------------------------------------------------------------------------------------
// §5. The simulated schedule.
// ---------------------------------------------------------------------------------------------

enum class JobKind : uint8_t { PRE, LEAF, FUSE, EXTRACT, COMBINE };

struct Job {
    JobKind kind{JobKind::LEAF};
    /// Tree node for `LEAF`/`FUSE`/`EXTRACT`; unused otherwise.
    uint32_t node{NONE};
    uint64_t priority{0};
    uint32_t pending{0};
    /// The dependency that finished last, i.e. the one that set `ready`.
    uint32_t ready_from{NONE};
    uint64_t ready{0};
    /// Bound core, or `-1` for "any". Fusions run on the core of their later-finishing child;
    /// extractions run on their root's core.
    int32_t core{-1};
    bool charge_dispatch{true};
    /// Filled by the run.
    uint64_t start{0};
    uint64_t finish{0};
    uint64_t ticks{0};
    /// `LEAF` only: the build half of `ticks`, so the critical path can charge `cp_build` and
    /// `cp_leaf` separately (amendment 1 §4). `ticks - build_ticks` is the solve half.
    uint64_t build_ticks{0};
    uint64_t dispatch_ticks{0};
    /// This job's own entry in `CellRunner::manager_chain`, when it was charged a dispatch.
    uint32_t manager_event{NONE};
    uint32_t predecessor{NONE};
    /// 0 = the dependency, 1 = the core, 2 = the manager.
    uint8_t determined_by{0};
    bool executed{false};
};

/// The whole of one shot's measurement, so that the row writer, the aggregate and the determinism
/// digest read the same values rather than three recomputations of them.
struct ShotMeasurement {
    uint32_t n_def{0};
    uint32_t n_edges{0};
    uint8_t escalated{0};
    /// 0 = none, 1 = a leaf truncated, 2 = a fusion truncated.
    uint8_t escalation_site{0};
    uint32_t n_pieces{0};
    uint32_t largest_piece{0};
    uint32_t n_refused{0};
    uint32_t n_fusions{0};
    uint32_t tree_depth{0};
    uint32_t max_fusion_boundary{0};
    uint32_t n_released{0};

    uint64_t ticks_uf{0};
    uint64_t ticks_tree{0};
    uint64_t ticks_scatter{0};
    uint64_t ticks_pre{0};
    uint64_t ticks_dispatch_total{0};
    /// The build half and the solve half of the leaves, separately: `_max` is the largest single
    /// leaf, `_sum` the total over all of them.
    uint64_t ticks_build_max{0};
    uint64_t ticks_build_sum{0};
    uint64_t ticks_solve_max{0};
    uint64_t ticks_solve_sum{0};
    uint64_t ticks_fuse_max{0};
    uint64_t ticks_fuse_sum{0};
    uint64_t ticks_extract_max{0};
    uint64_t ticks_combine{0};

    uint64_t cp_manager{0};
    uint64_t cp_build{0};
    uint64_t cp_leaf{0};
    uint64_t cp_fuse{0};
    uint64_t cp_extract{0};
    uint64_t cp_idle{0};
    uint32_t cp_len{0};

    uint64_t sparse_k{0};
    uint64_t fallback{0};
    uint64_t system{0};

    /// The fused answer, for `--verify`. Not a column.
    pm::obs_int obs_mask{0};
    pm::total_weight_int weight{0};
};

struct Row {
    uint64_t shot;
    ShotMeasurement m;
};

struct CellLabels {
    size_t distance{0};
    double noise{0};
    double horizon{0};
    size_t cores{0};
    double alpha{0};
};

struct RowArena {
    std::vector<Row> rows;
    CellLabels labels;
    std::ofstream* out{nullptr};
    uint64_t rows_this_cell{0};
    uint64_t rows_total{0};
    int64_t truncated_at_shot{-1};

    void configure() {
        rows.reserve(ROW_ARENA_CAPACITY);
    }

    inline void append(const Row& row) {
        if (rows.size() >= ROW_ARENA_CAPACITY) {
            if (truncated_at_shot < 0)
                truncated_at_shot = (int64_t)row.shot;
            return;
        }
        rows.push_back(row);
    }

    /// Called from between-shot positions only, so no timed region is ever open across it.
    void flush() {
        for (const Row& row : rows) {
            const ShotMeasurement& m = row.m;
            *out << row.shot << "," << m.n_def << "," << m.n_edges << "," << (int)m.escalated << "," << m.n_pieces
                 << "," << m.largest_piece << "," << m.n_refused << "," << m.n_fusions << "," << m.tree_depth << ","
                 << m.max_fusion_boundary << "," << m.n_released << "," << m.ticks_uf << "," << m.ticks_tree << ","
                 << m.ticks_scatter << "," << m.ticks_pre << "," << m.ticks_dispatch_total << "," << m.ticks_build_max
                 << "," << m.ticks_build_sum << "," << m.ticks_solve_max << "," << m.ticks_solve_sum << ","
                 << m.ticks_fuse_max << "," << m.ticks_fuse_sum << "," << m.ticks_extract_max << "," << m.ticks_combine
                 << "," << m.cp_manager << "," << m.cp_build << "," << m.cp_leaf << "," << m.cp_fuse << ","
                 << m.cp_extract << "," << m.cp_idle << "," << m.cp_len << "," << m.sparse_k << "," << m.fallback << ","
                 << m.system << "\n";
        }
        rows_this_cell += rows.size();
        rows_total += rows.size();
        rows.clear();
        out->flush();
    }
};

struct TickSeries {
    __int128 sum{0};
    __int128 sum_sq{0};
    uint64_t max{0};
    std::vector<long long> samples;

    void reserve(size_t count) {
        samples.reserve(count);
    }
    void clear() {
        sum = 0;
        sum_sq = 0;
        max = 0;
        samples.clear();
    }
    void add(uint64_t ticks) {
        sum += (__int128)ticks;
        sum_sq += (__int128)ticks * (__int128)ticks;
        max = std::max(max, ticks);
        samples.push_back((long long)ticks);
    }
    double mean() const {
        return samples.empty() ? 0.0 : (double)sum / (double)samples.size();
    }
};

/// One `(d, p, T, k)` cell of `agg.json`.
struct CellAggregate {
    size_t distance{0};
    double noise{0};
    double horizon{0};
    size_t cores{0};
    double alpha{0};
    uint64_t shots{0};
    uint64_t escalations{0};
    uint64_t escalations_in_leaf{0};
    uint64_t escalations_in_fusion{0};

    TickSeries fallback;
    TickSeries sparse_k;
    TickSeries system;

    /// Means of the `ticks_*` and `cp_*` columns, over the **non-escalating** shots (deviation 5).
    uint64_t module_shots{0};
    __int128 uf_sum{0};
    __int128 tree_sum{0};
    __int128 scatter_sum{0};
    __int128 pre_sum{0};
    __int128 dispatch_sum{0};
    __int128 build_max_sum{0};
    __int128 build_sum_sum{0};
    __int128 solve_max_sum{0};
    __int128 solve_sum_sum{0};
    __int128 fuse_max_sum{0};
    __int128 fuse_sum_sum{0};
    __int128 extract_max_sum{0};
    __int128 combine_sum{0};
    __int128 cp_manager_sum{0};
    __int128 cp_build_sum{0};
    __int128 cp_leaf_sum{0};
    __int128 cp_fuse_sum{0};
    __int128 cp_extract_sum{0};
    __int128 cp_idle_sum{0};
    __int128 cp_len_sum{0};

    std::vector<uint64_t> n_pieces_hist;
    std::vector<uint64_t> tree_depth_hist;
    std::vector<uint64_t> max_fusion_boundary_hist;
    std::vector<uint64_t> n_released_hist;

    int64_t raw_truncated_at_shot{-1};
    /// `1` / `0` / `-1` for pass / fail / not run.
    int determinism_check{-1};
    /// §11 tests 3 and 4, checked on every measured shot when asserts are live.
    int64_t threshold_violations{0};
    int64_t cut_on_heaviest_violations{0};
    int64_t feasibility_violations{0};
    /// Amendment 1 tests 7, 8 and 9, debug builds only.
    int64_t manager_graph_write_violations{0};
    int64_t build_independence_violations{0};
    int64_t built_before_fused_violations{0};
    int64_t verify_shots{0};
    /// §8's escalation clause, counted rather than asserted; see deviation 9.
    int64_t escalation_divergence_fused_only{0};
    int64_t escalation_divergence_fused_only_in_leaf{0};
    int64_t escalation_divergence_fused_only_in_fusion{0};
    int64_t escalation_divergence_monolithic_only{0};

    void reset(size_t histogram_size) {
        shots = 0;
        escalations = 0;
        escalations_in_leaf = 0;
        escalations_in_fusion = 0;
        fallback.clear();
        sparse_k.clear();
        system.clear();
        module_shots = 0;
        uf_sum = tree_sum = scatter_sum = pre_sum = dispatch_sum = 0;
        build_max_sum = build_sum_sum = solve_max_sum = solve_sum_sum = 0;
        fuse_max_sum = fuse_sum_sum = extract_max_sum = combine_sum = 0;
        cp_manager_sum = cp_build_sum = cp_leaf_sum = cp_fuse_sum = cp_extract_sum = cp_idle_sum = cp_len_sum = 0;
        n_pieces_hist.assign(histogram_size, 0);
        tree_depth_hist.assign(histogram_size, 0);
        max_fusion_boundary_hist.assign(histogram_size, 0);
        n_released_hist.assign(histogram_size, 0);
        raw_truncated_at_shot = -1;
        determinism_check = -1;
        threshold_violations = 0;
        cut_on_heaviest_violations = 0;
        feasibility_violations = 0;
        manager_graph_write_violations = 0;
        build_independence_violations = 0;
        built_before_fused_violations = 0;
        verify_shots = 0;
        escalation_divergence_fused_only = 0;
        escalation_divergence_fused_only_in_leaf = 0;
        escalation_divergence_fused_only_in_fusion = 0;
        escalation_divergence_monolithic_only = 0;
    }
};

bool stderr_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    if (isatty(fileno(stderr)) != 1)
        return false;
    const char* term = std::getenv("TERM");
    return term != nullptr && std::string(term) != "dumb";
#else
    return false;
#endif
}

constexpr size_t PROGRESS_BAR_WIDTH = 24;

void draw_progress(
    bool interactive, const char* label, size_t cell_index, size_t total_cells, size_t done, size_t cell_shots) {
    double fraction = cell_shots ? (double)done / (double)cell_shots : 1.0;
    if (!interactive) {
        std::fprintf(
            stderr, "[%zu/%zu] %s  %5.1f%%  %zu/%zu shots\n", cell_index, total_cells, label, 100.0 * fraction, done,
            cell_shots);
        std::fflush(stderr);
        return;
    }
    char bar[PROGRESS_BAR_WIDTH + 1];
    size_t filled = (size_t)(fraction * (double)PROGRESS_BAR_WIDTH);
    for (size_t i = 0; i < PROGRESS_BAR_WIDTH; i++)
        bar[i] = i < filled ? '#' : '.';
    bar[PROGRESS_BAR_WIDTH] = '\0';
    std::fprintf(
        stderr, "\r[%zu/%zu] %s  [%s] %5.1f%%  %zu/%zu shots\x1b[K", cell_index, total_cells, label, bar,
        100.0 * fraction, done, cell_shots);
    std::fflush(stderr);
}

void draw_setup(bool interactive, const char* label, size_t cell_index, size_t total_cells) {
    std::fprintf(
        stderr,
        interactive ? "\r[%zu/%zu] %s  building DEM, ball tables and the shared instance...\x1b[K"
                    : "[%zu/%zu] %s  building DEM, ball tables and the shared instance...\n",
        cell_index, total_cells, label);
    std::fflush(stderr);
}

/// The generator call, as one string, so `run.log` records what `(d, p)` actually meant. It is
/// `profiler_util.h`'s `ShotSampler::make` — the same call the rest of this corpus is generated by,
/// with the two deviations that call carries (see the header).
std::string generator_call(size_t distance, size_t rounds, double noise) {
    std::ostringstream out;
    out << "stim::CircuitGenParameters(rounds=" << rounds << ", distance=" << distance
        << ", task=\"rotated_memory_x\"); after_clifford_depolarization=" << noise
        << "; after_reset_flip_probability=" << noise << "; before_measure_flip_probability=" << noise
        << "; before_round_data_depolarization=0 (unset)"
        << "; stim::generate_surface_code_circuit(...).circuit"
        << "; stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false)";
    return out.str();
}

std::string stamp() {
    std::time_t now = std::time(nullptr);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
    return std::string(buffer);
}

/// The whole of one `(d, p, T, k)` cell's per-shot work, and the only place a clock is read.
struct CellRunner {
    BallDecoder* decoder{nullptr};
    pm::Mwpm* fallback_instance{nullptr};
    Cut* cut{nullptr};
    /// The shared instance: one `BallMwpm` holding the shot's whole `H`.
    BallMwpm* shared{nullptr};
    fus::FusionInstance* fusion{nullptr};
    /// §8's monolithic reference, on its own instance so the two solves share no state.
    BallMwpm* reference{nullptr};
    size_t k{0};

    /// The adjacency index of each crossing edge's two directions. Each half is filled by the build
    /// of the piece that owns that endpoint (amendment 1 §3.2), where the index is known, rather
    /// than binary-searched for on the manager core.
    std::vector<uint32_t> cross_index_lo;
    std::vector<uint32_t> cross_index_hi;
    /// The observable mask of every edge of `H`, derived from the ball tables once per shot with the
    /// rest of the input and outside every timed region. The design's input record is
    /// `(u, v, w, obs_mask)`; this is the `obs_mask` column of it.
    std::vector<pm::obs_int> edge_obs;
    /// One byte per defect: whether its adjacency has been written by its piece's build. Cleared
    /// with the node pool between shots, outside every timed region.
    std::vector<uint8_t> built;
    /// The solver core's gather buffers for `build_piece`.
    fus::PieceBuildScratch build_scratch;
    /// The build's read-only view of the shot, rebound once per shot outside every timed region.
    fus::PieceBuild build_input;

    /// §5's job list and core state.
    std::vector<Job> jobs;
    std::vector<uint32_t> ready_heap;
    std::vector<uint64_t> core_free;
    std::vector<uint32_t> core_last_job;
    /// Job id per tree node, so a finished child can find its parent's job.
    std::vector<uint32_t> job_of_node;
    std::vector<uint32_t> extract_of_node;
    /// Per-root extraction results, reduced by `COMBINE`.
    std::vector<pm::MatchingResult> root_results;

    /// The fallback's output buffer. Zeroed before the region, never inside it.
    std::vector<uint8_t> obs_bytes;
    /// All of `H`'s defects, for the monolithic reference.
    std::vector<uint32_t> all_defects;

    uint64_t manager_free{0};
    uint32_t manager_last_job{NONE};
    /// The manager core's own timeline, as a contiguous chain of its jobs: entry 0 is `PRE` and
    /// every later entry is one charged dispatch (and, last, `COMBINE` if the manager was what held
    /// it up). It exists because a manager-determined start cannot chain back to the *job* that was
    /// dispatched last — that job is still running on a solver core, and its finish is in the
    /// future — but it can chain back to the manager's own previous piece of work, which is exactly
    /// what was in the way.
    std::vector<uint64_t> manager_chain;

    void configure(size_t n_max, size_t num_observables, size_t k_max) {
        obs_bytes.assign(std::max<size_t>(1, num_observables), 0);
        core_free.assign(k_max, 0);
        core_last_job.assign(k_max, NONE);
        edge_obs.reserve(INITIAL_EDGE_CAPACITY);
        built.assign(n_max, 0);
        build_scratch.reserve(n_max, 4 * INITIAL_EDGE_CAPACITY);
        cross_index_lo.reserve(INITIAL_EDGE_CAPACITY);
        cross_index_hi.reserve(INITIAL_EDGE_CAPACITY);
        jobs.reserve(4 * n_max + 8);
        ready_heap.reserve(4 * n_max + 8);
        job_of_node.assign(2 * n_max + 1, NONE);
        extract_of_node.assign(2 * n_max + 1, NONE);
        root_results.reserve(n_max);
        all_defects.reserve(n_max);
        walk_stack.reserve(2 * n_max + 2);
        manager_chain.reserve(4 * n_max + 8);
    }

    /// Total capacity across every buffer a measured shot could make grow, split three ways so that
    /// a `no_growth_after_warmup=0` says *what* grew rather than only that something did.
    ///
    ///  - `input` is the ball-table `H` builder's own arena, which this binary shares with the rest
    ///    of the corpus and does not own;
    ///  - `solver` is the shared instance's node pool;
    ///  - `own` is everything §3 and §5 allocate: the cut, the mask, the job list and the scatter.
    struct CapacityBreakdown {
        uint64_t input{0};
        uint64_t solver{0};
        uint64_t own{0};

        uint64_t total() const {
            return input + solver + own;
        }
        bool operator==(const CapacityBreakdown& other) const {
            return input == other.input && solver == other.solver && own == other.own;
        }
    };

    CapacityBreakdown capacity_breakdown() const {
        CapacityBreakdown out;
        out.input = decoder->arena.grow_events;
        out.solver = shared->grow_events;
        out.own = (uint64_t)edge_obs.capacity() + built.capacity() + build_scratch.capacity_signature() +
                  cross_index_lo.capacity() + cross_index_hi.capacity() + jobs.capacity() + ready_heap.capacity() +
                  root_results.capacity() + all_defects.capacity() + walk_stack.capacity() +
                  manager_chain.capacity() + refused_slot.capacity() + fusion->mask.capacity_signature() +
                  cut->capacity_signature();
        return out;
    }

    // -----------------------------------------------------------------------------------------
    // §3. Preprocessing, on the manager core.
    // -----------------------------------------------------------------------------------------

    /// Amendment 1 §1. The hardware's sort, stood in for by a `std::sort`, **outside every timed
    /// region**: the edge list reaches the solver already in non-decreasing weight, ties by
    /// `(u, v)`, because the sort happens in the same place the edge list is produced.
    ///
    /// Leaves `cut->order` a permutation of the edge indices in that order. `h.edges` itself is left
    /// alone — it is sorted by `(i, j)`, which is what `--verify`'s monolithic `rebuild` needs.
    void sort_edges(const BallGraph& h) {
        Cut& c = *cut;
        c.order.resize(h.edges.size());
        for (uint32_t e = 0; e < (uint32_t)h.edges.size(); e++)
            c.order[e] = e;
        // `(w, u, v)` is a strict total order on `H`'s edges, so there is no tie for an unstable
        // sort to break differently from one run to the next.
        std::sort(c.order.begin(), c.order.end(), [&](uint32_t a, uint32_t b) {
            const BallGraphEdge& x = h.edges[a];
            const BallGraphEdge& y = h.edges[b];
            if (x.w_int != y.w_int)
                return x.w_int < y.w_int;
            if (x.i != y.i)
                return x.i < y.i;
            return x.j < y.j;
        });
    }

    /// §3.2. Bounded union-find with the boundary weights, then (deviation 3) the quotient graph
    /// from the pieces the pass settled on.
    void union_find(const BallGraph& h) {
        Cut& c = *cut;
        uint32_t n = (uint32_t)h.num_nodes();
        uint32_t threshold = (uint32_t)((n + k - 1) / k);
        if (threshold == 0)
            threshold = 1;

        for (uint32_t i = 0; i < n; i++) {
            c.parent[i] = i;
            c.set_size[i] = 1;
            c.w_boundary[i] = WEIGHT_INF;
            c.bnd[i] = WEIGHT_INF;
            c.via_dummy[i] = 0;
        }
        // `H` gives a defect at most one boundary edge, so this is a fill rather than a reduction.
        for (const BallBoundaryEdge& edge : h.boundary_edges) {
            c.w_boundary[edge.i] = edge.w_int;
            c.bnd[edge.i] = edge.w_int;
        }

        c.refused.clear();
        for (uint32_t position = 0; position < (uint32_t)c.order.size(); position++) {
            uint32_t e = c.order[position];
            const BallGraphEdge& edge = h.edges[e];
            uint32_t a = c.find(edge.i);
            uint32_t b = c.find(edge.j);
            if (a == b)
                continue;
            if ((uint64_t)c.set_size[a] + c.set_size[b] <= threshold) {
                // Weighted union, ties to the smaller index, so the link is a function of `H` alone.
                if (c.set_size[a] < c.set_size[b] || (c.set_size[a] == c.set_size[b] && b < a))
                    std::swap(a, b);
                c.parent[b] = a;
                c.set_size[a] += c.set_size[b];
                continue;
            }
            // Refused: this edge is a crossing edge, and both endpoints gain (or tighten) a dummy
            // boundary at half its weight. Edges arrive in non-decreasing weight and `floor(w / 2)`
            // is monotone in `w`, so the first refusal at a defect already sets its final `bnd` and
            // the two lines below are no-ops on every later one.
            c.refused.push_back(e);
            pm::weight_int half = fus::even_dummy_weight(edge.w_int);
            if (half < c.bnd[edge.i]) {
                c.bnd[edge.i] = half;
                c.via_dummy[edge.i] = 1;
            }
            if (half < c.bnd[edge.j]) {
                c.bnd[edge.j] = half;
                c.via_dummy[edge.j] = 1;
            }
        }

        // Pieces, numbered in ascending minimum member. Iterating `v` upwards, the first `v` of a
        // set *is* its minimum, so the numbering falls out of the pass with no sort.
        c.epoch++;
        c.n_pieces = 0;
        c.largest_piece = 0;
        for (uint32_t v = 0; v < n; v++) {
            uint32_t root = c.find(v);
            if (c.root_stamp[root] != c.epoch) {
                c.root_stamp[root] = c.epoch;
                c.piece_index_of_root[root] = c.n_pieces;
                c.n_pieces++;
            }
            c.piece_of[v] = c.piece_index_of_root[root];
        }
        // Member blocks, ascending within each piece.
        for (uint32_t i = 0; i <= c.n_pieces; i++)
            c.piece_offsets[i] = 0;
        for (uint32_t v = 0; v < n; v++)
            c.piece_offsets[c.piece_of[v] + 1]++;
        for (uint32_t i = 0; i < c.n_pieces; i++) {
            c.piece_offsets[i + 1] += c.piece_offsets[i];
            c.fill_cursor[i] = c.piece_offsets[i];
        }
        for (uint32_t v = 0; v < n; v++)
            c.piece_members[c.fill_cursor[c.piece_of[v]]++] = v;
        for (uint32_t i = 0; i < c.n_pieces; i++)
            c.largest_piece = std::max(c.largest_piece, c.piece_offsets[i + 1] - c.piece_offsets[i]);

        // Deviation 3: the quotient graph, from the refused edges' final piece ids.
        c.size_quotient_table(c.refused.size());
        c.qslots.clear();
        for (uint32_t e : c.refused) {
            const BallGraphEdge& edge = h.edges[e];
            uint32_t pa = c.piece_of[edge.i];
            uint32_t pb = c.piece_of[edge.j];
            uint32_t lo = std::min(pa, pb);
            uint32_t hi = std::max(pa, pb);
            uint64_t key = ((uint64_t)lo << 32) | hi;
            size_t slot = (size_t)((key * 1099511628211ull) >> 24) & (c.qcap - 1);
            while (c.qkey[slot] != UINT64_MAX && c.qkey[slot] != key)
                slot = (slot + 1) & (c.qcap - 1);
            if (c.qkey[slot] == UINT64_MAX) {
                c.qkey[slot] = key;
                c.qcount[slot] = 0;
                c.qslots.push_back((uint32_t)slot);
            }
            c.qcount[slot]++;
        }
    }

    /// §3.3. The cheapest-boundary-first fusion tree, then (deviation 4) each fusion's interior
    /// crossing edges.
    void build_tree(const BallGraph& h) {
        Cut& c = *cut;
        c.n_tree_nodes = c.n_pieces;
        c.n_fusions = 0;
        c.max_fusion_boundary = 0;
        c.depth = 0;
        for (uint32_t i = 0; i < c.n_pieces; i++) {
            c.tree_child0[i] = NONE;
            c.tree_child1[i] = NONE;
            c.tree_parent[i] = NONE;
            c.tree_cost[i] = 0;
            c.tree_alive[i] = 1;
            c.tuf[i] = i;
        }
        for (uint32_t node : c.adj_touched)
            c.adj[node].clear();
        c.adj_touched.clear();
        c.heap.clear();

        for (uint32_t slot : c.qslots) {
            uint64_t key = c.qkey[slot];
            uint32_t lo = (uint32_t)(key >> 32);
            uint32_t hi = (uint32_t)(key & 0xffffffffu);
            uint32_t count = c.qcount[slot];
            if (c.adj[lo].empty())
                c.adj_touched.push_back(lo);
            if (c.adj[hi].empty())
                c.adj_touched.push_back(hi);
            c.adj[lo].push_back(QuotientAdj{hi, count});
            c.adj[hi].push_back(QuotientAdj{lo, count});
            c.heap.push_back(MergeCandidate{count, lo, hi});
            c.qkey[slot] = UINT64_MAX;  // The table is left clean for the next shot.
        }
        c.qslots.clear();
        std::make_heap(c.heap.begin(), c.heap.end());

        while (!c.heap.empty()) {
            std::pop_heap(c.heap.begin(), c.heap.end());
            MergeCandidate candidate = c.heap.back();
            c.heap.pop_back();
            // Stale: an endpoint has already been merged into something.
            if (!c.tree_alive[candidate.lo] || !c.tree_alive[candidate.hi])
                continue;

            uint32_t a = candidate.lo;
            uint32_t b = candidate.hi;
            uint32_t p = c.n_tree_nodes++;
            c.tree_child0[p] = a;
            c.tree_child1[p] = b;
            c.tree_parent[p] = NONE;
            c.tree_parent[a] = p;
            c.tree_parent[b] = p;
            c.tree_cost[p] = candidate.count;
            c.tree_alive[a] = 0;
            c.tree_alive[b] = 0;
            c.tree_alive[p] = 1;
            c.tuf[p] = p;
            c.tuf[a] = p;
            c.tuf[b] = p;
            c.n_fusions++;
            c.max_fusion_boundary = std::max(c.max_fusion_boundary, candidate.count);

            // Small-to-large: the larger child's list becomes the parent's, and the smaller one is
            // folded in. Entries naming a node that has since been merged are resolved through
            // `tuf`; entries resolving to the parent itself are the `a`-`b` edge and are dropped.
            std::vector<QuotientAdj>& big = c.adj[a].size() >= c.adj[b].size() ? c.adj[a] : c.adj[b];
            std::vector<QuotientAdj>& small = c.adj[a].size() >= c.adj[b].size() ? c.adj[b] : c.adj[a];
            if (c.adj[p].empty())
                c.adj_touched.push_back(p);
            c.adj[p].swap(big);
            c.mark_epoch++;
            size_t write = 0;
            for (size_t i = 0; i < c.adj[p].size(); i++) {
                uint32_t other = c.find_tree(c.adj[p][i].other);
                if (other == p)
                    continue;
                if (c.mark_stamp[other] == c.mark_epoch) {
                    c.adj[p][c.mark_pos[other]].count += c.adj[p][i].count;
                    continue;
                }
                c.mark_stamp[other] = c.mark_epoch;
                c.mark_pos[other] = (uint32_t)write;
                c.adj[p][write++] = QuotientAdj{other, c.adj[p][i].count};
            }
            c.adj[p].resize(write);
            for (const QuotientAdj& entry : small) {
                uint32_t other = c.find_tree(entry.other);
                if (other == p)
                    continue;
                if (c.mark_stamp[other] == c.mark_epoch) {
                    c.adj[p][c.mark_pos[other]].count += entry.count;
                    continue;
                }
                c.mark_stamp[other] = c.mark_epoch;
                c.mark_pos[other] = (uint32_t)c.adj[p].size();
                c.adj[p].push_back(QuotientAdj{other, entry.count});
            }
            big.clear();
            small.clear();
            for (const QuotientAdj& entry : c.adj[p]) {
                c.heap.push_back(
                    MergeCandidate{entry.count, std::min(p, entry.other), std::max(p, entry.other)});
                std::push_heap(c.heap.begin(), c.heap.end());
            }
        }

        // Roots, ascending, and the tree's depth. Depths are computed top-down over the internal
        // nodes, which were created in increasing id order, so a parent is always seen first when
        // the list is walked backwards.
        c.tree_roots.clear();
        for (uint32_t node = 0; node < c.n_tree_nodes; node++) {
            c.tree_depth[node] = 0;
            if (c.tree_alive[node])
                c.tree_roots.push_back(node);
        }
        for (uint32_t node = c.n_tree_nodes; node-- > c.n_pieces;) {
            uint32_t depth = c.tree_depth[node] + 1;
            c.tree_depth[c.tree_child0[node]] = depth;
            c.tree_depth[c.tree_child1[node]] = depth;
            c.depth = std::max(c.depth, depth);
        }

        // Deviation 4: each refused edge becomes interior at the lowest common ancestor of its two
        // pieces' leaves, which is by definition the fusion that unmasks it.
        for (uint32_t i = 0; i <= c.n_tree_nodes; i++)
            c.interior_offsets[i] = 0;
        c.interior_edges.resize(c.refused.size());
        // Two pieces joined by a refused edge are joined by a quotient edge, and the greedy does not
        // stop until the heap is empty, so they always end up under one root and the walk always
        // meets. The `NONE` guards are there so that a bug is a thrown error rather than a read off
        // the end of `tree_parent`.
        auto lca = [&](uint32_t x, uint32_t y) -> uint32_t {
            while (x != NONE && y != NONE && c.tree_depth[x] > c.tree_depth[y])
                x = c.tree_parent[x];
            while (x != NONE && y != NONE && c.tree_depth[y] > c.tree_depth[x])
                y = c.tree_parent[y];
            while (x != NONE && y != NONE && x != y) {
                x = c.tree_parent[x];
                y = c.tree_parent[y];
            }
            if (x == NONE || y == NONE || x != y)
                throw std::runtime_error("fusion tree: a crossing edge's pieces have no common ancestor");
            return x;
        };
        for (uint32_t e : c.refused) {
            uint32_t node = lca(c.piece_of[h.edges[e].i], c.piece_of[h.edges[e].j]);
            c.interior_offsets[node + 1]++;
        }
        for (uint32_t i = 0; i < c.n_tree_nodes; i++) {
            c.interior_offsets[i + 1] += c.interior_offsets[i];
            c.fill_cursor[i] = c.interior_offsets[i];
        }
        for (uint32_t e : c.refused) {
            uint32_t node = lca(c.piece_of[h.edges[e].i], c.piece_of[h.edges[e].j]);
            c.interior_edges[c.fill_cursor[node]++] = e;
        }
    }

    /// Amendment 1 §2.2. A light pass, and no solver graph work at all: the sorted edge list grouped
    /// by piece so that each leaf can build its own part of the graph without scanning the whole
    /// list, plus the bookkeeping a fusion reads.
    ///
    /// It reads no ball table and touches no solver node. It is `O(n + m)` integer work in fixed
    /// arenas, and amendment 1 test 7 checks the "no solver node" half of that outright.
    void scatter(const BallGraph& h) {
        Cut& c = *cut;
        uint32_t n = (uint32_t)h.num_nodes();
        uint32_t m = (uint32_t)h.edges.size();

        // Which entry of the crossing-edge list each refused edge is. Inside this timer rather than
        // in the gap between `tree` and `scatter` where it used to sit untimed; it is O(refused).
        refused_slot.resize(std::max<size_t>(1, m));
        for (uint32_t r = 0; r < (uint32_t)c.refused.size(); r++)
            refused_slot[c.refused[r]] = r;

        // ---- Step 1: edge -> piece, and each defect's adjacency degree.
        //
        // The piece of an edge is its endpoints' final union-find root, which is why this runs after
        // `uf`. An interior edge belongs to the one piece; a crossing edge belongs to **both**
        // endpoint pieces, once each, since each endpoint writes its own adjacency entry.
        for (uint32_t i = 0; i <= c.n_pieces; i++)
            c.piece_edge_offsets[i] = 0;
        // The boundary half-edge, where there is one, is adjacency entry 0. `bnd` only ever rises at
        // a fusion, and only as far as `w_boundary`, so a defect that has no slot here never gains
        // one and a defect that has one never loses it.
        for (uint32_t v = 0; v < n; v++)
            c.node_degree[v] = c.bnd[v] == WEIGHT_INF ? 0 : 1;
        for (uint32_t e = 0; e < m; e++) {
            const BallGraphEdge& edge = h.edges[e];
            uint32_t pa = c.piece_of[edge.i];
            uint32_t pb = c.piece_of[edge.j];
            c.piece_edge_offsets[pa + 1]++;
            if (pb != pa)
                c.piece_edge_offsets[pb + 1]++;
            c.node_degree[edge.i]++;
            c.node_degree[edge.j]++;
        }

        // ---- Step 2: counting sort by piece. The fill walks the sorted order, so each piece's
        // slice comes out in non-decreasing weight, exactly as the whole list is.
        for (uint32_t i = 0; i < c.n_pieces; i++) {
            c.piece_edge_offsets[i + 1] += c.piece_edge_offsets[i];
            c.fill_cursor[i] = c.piece_edge_offsets[i];
        }
        c.piece_edges.resize(c.piece_edge_offsets[c.n_pieces]);
        for (uint32_t position = 0; position < m; position++) {
            uint32_t e = c.order[position];
            const BallGraphEdge& edge = h.edges[e];
            uint32_t pa = c.piece_of[edge.i];
            uint32_t pb = c.piece_of[edge.j];
            c.piece_edges[c.fill_cursor[pa]++] = e;
            if (pb != pa)
                c.piece_edges[c.fill_cursor[pb]++] = e;
        }

        // ---- Step 3: the per-piece defect blocks are filled by `uf`, where the pass that numbers
        // the pieces already has them in hand.

        // ---- Step 4: the bookkeeping a fusion reads.
        //
        // The crossing edges incident to each defect, so that a fusion can recompute `wmin` over the
        // ones that are still external without walking the whole refused list.
        for (uint32_t i = 0; i <= n; i++)
            c.cross_offsets[i] = 0;
        for (uint32_t e : c.refused) {
            c.cross_offsets[h.edges[e].i + 1]++;
            c.cross_offsets[h.edges[e].j + 1]++;
        }
        for (uint32_t i = 0; i < n; i++) {
            c.cross_offsets[i + 1] += c.cross_offsets[i];
            c.fill_cursor[i] = c.cross_offsets[i];
        }
        c.cross_edges.resize(c.cross_offsets[n]);
        for (uint32_t e : c.refused) {
            c.cross_edges[c.fill_cursor[h.edges[e].i]++] = e;
            c.cross_edges[c.fill_cursor[h.edges[e].j]++] = e;
        }
        cross_index_lo.resize(c.refused.size());
        cross_index_hi.resize(c.refused.size());

        // The mask arena's offsets, which the builds write their mask bytes through and which a
        // fusion unmasks through. They are a prefix sum of the degrees counted above rather than of
        // the instance's adjacency, because no adjacency exists yet — that is the point of the
        // amendment. The build scratch shares the same offsets, so two pieces can never gather into
        // the same word.
        fusion->size_mask_for_shot(c.node_degree.data(), n);
        build_scratch.size_for_shot(fusion->mask.total_entries(), n);

        // The runtime replay of the tree's unions starts from the leaves again.
        for (uint32_t node = 0; node < c.n_tree_nodes; node++)
            c.tuf[node] = node;
    }

    /// Amendment 1 §3.2, on this core rather than the manager's: one piece's defects, and only
    /// those, get their adjacency, weights, observable masks, mask bytes and boundary edge.
    void build_piece(uint32_t piece) {
        Cut& c = *cut;
        fus::build_piece(
            *fusion, build_scratch, build_input, piece, c.piece_members.data() + c.piece_offsets[piece],
            c.piece_offsets[piece + 1] - c.piece_offsets[piece],
            c.piece_edges.data() + c.piece_edge_offsets[piece],
            c.piece_edge_offsets[piece + 1] - c.piece_edge_offsets[piece]);
    }

    // -----------------------------------------------------------------------------------------
    // §4.4. One fusion.
    // -----------------------------------------------------------------------------------------

    /// Fuses tree node `p`'s two children. Returns the status and adds to `released`.
    TimelineStatus fuse(const BallGraph& h, uint32_t p, uint32_t& released) {
        Cut& c = *cut;
        check_built_before_fused(p);
        // A fusion is the other place growth starts, so the clock has to be even before a region is
        // released into it. See `fusion::align_clock_for_new_growth`.
        fus::align_clock_for_new_growth(shared->mwpm);
        // Union first, so that "is `v` inside `p`?" is a `find_tree` away for the rest of this.
        c.tuf[c.tree_child0[p]] = p;
        c.tuf[c.tree_child1[p]] = p;

        uint32_t begin = c.interior_offsets[p];
        uint32_t end = c.interior_offsets[p + 1];

        // Step 1: every crossing edge of `p` becomes a real edge again.
        for (uint32_t at = begin; at < end; at++) {
            uint32_t e = c.interior_edges[at];
            size_t r = refused_position(e);
            const BallGraphEdge& edge = h.edges[e];
            fus::unmask_crossing_edge(*fusion, edge.i, cross_index_lo[r], edge.j, cross_index_hi[r]);
        }

        // Step 2: the dummies those edges were standing in for are raised, and every region a raised
        // dummy was holding matched is released.
        c.defect_epoch++;
        for (uint32_t at = begin; at < end; at++) {
            const BallGraphEdge& edge = h.edges[c.interior_edges[at]];
            for (uint32_t u : {edge.i, edge.j}) {
                if (c.defect_stamp[u] == c.defect_epoch || !c.via_dummy[u])
                    continue;
                c.defect_stamp[u] = c.defect_epoch;
                // `wmin` over the crossing edges at `u` that still leave `p`.
                pm::weight_int wmin = WEIGHT_INF;
                for (uint32_t at2 = c.cross_offsets[u]; at2 < c.cross_offsets[u + 1]; at2++) {
                    const BallGraphEdge& other = h.edges[c.cross_edges[at2]];
                    uint32_t far = other.i == u ? other.j : other.i;
                    if (c.find_tree(c.piece_of[far]) == p)
                        continue;
                    wmin = std::min(wmin, other.w_int);
                }
                pm::weight_int half = wmin == WEIGHT_INF ? WEIGHT_INF : fus::even_dummy_weight(wmin);
                pm::weight_int want = std::min(c.w_boundary[u], half);
                bool dummy = want < c.w_boundary[u];
                // The bookkeeping is updated whether or not the instance is touched: `want` can
                // equal the current `bnd` and still stop being a dummy, when the last crossing edge
                // goes and `floor(wmin/2)` happened to equal the real boundary's own weight. Leaving
                // `via_dummy` set there would keep an observable mask of 0 on an edge that is now
                // the code boundary.
                pm::weight_int previous = c.bnd[u];
                bool was_dummy = c.via_dummy[u] != 0;
                c.bnd[u] = want;
                c.via_dummy[u] = dummy ? 1 : 0;
                if (want == previous && dummy == was_dummy)
                    continue;
                pm::weight_int weight = want == WEIGHT_INF ? fusion->unreachable_boundary_weight() : want;
                // A dummy carries no observable, and neither does an unreachable edge; only a real
                // code boundary does.
                pm::obs_int obs = (dummy || want == WEIGHT_INF) ? (pm::obs_int)0 : c.obs_boundary[u];
                if (fus::set_boundary_and_release(shared->mwpm, u, weight, obs))
                    released++;
            }
        }

        // Step 3.
        return fus::run_until_settled(shared->mwpm);
    }

    /// The index of edge `e` within `cut->refused`. The refused list is in the input's sorted weight
    /// order and the interior lists are built from it, so this is a lookup into a small side table
    /// rather than a search: `refused_slot[e]` is filled at scatter time, and a build reads it to
    /// record where it put a crossing edge's adjacency entry.
    std::vector<uint32_t> refused_slot;
    inline size_t refused_position(uint32_t e) const {
        return refused_slot[e];
    }

    // -----------------------------------------------------------------------------------------
    // §5. One shot.
    // -----------------------------------------------------------------------------------------

    template <bool Profile>
    void run_shot(const std::vector<uint64_t>& shot, ShotMeasurement& out) {
        Cut& c = *cut;
        auto tick = []() -> uint64_t {
            if constexpr (Profile)
                return hires_now_ns();
            else
                return 0;
        };

        // ---- The input. Producing the edge list is outside every timed region.
        decoder->compute_seeded_detection_events(shot, decoder->seeded_scratch);
        build_ball_graph(
            decoder->tables, decoder->seeded_scratch, decoder->horizon, decoder->arena, decoder->config.mode, nullptr,
            nullptr);
        const BallGraph& h = decoder->arena.graph;
        uint32_t n = (uint32_t)h.num_nodes();

        out = ShotMeasurement();
        out.n_def = n;
        out.n_edges = (uint32_t)h.edges.size();

        // `obs_boundary` is read inside `fuse`, and `edge_obs` by every build, so both are filled
        // here, with the input, rather than re-derived from the ball tables inside a timed region.
        // The design's input record is `(u, v, w, obs_mask)`; this is where the mask column of it
        // comes from.
        for (const BallBoundaryEdge& edge : h.boundary_edges)
            c.obs_boundary[edge.i] = boundary_obs_mask(edge.det);
        edge_obs.resize(h.edges.size());
        for (size_t e = 0; e < h.edges.size(); e++)
            edge_obs[e] = edge_obs_mask(h.edges[e].entry);

        // The hardware's sort, outside every timed region (amendment 1 §1).
        sort_edges(h);

        // The previous shot's adjacency is what makes a node "built"; dropping it is what marks
        // every node unbuilt again. Outside every timed region, and never charged to any core.
        // `H`'s defects are always `0 .. n - 1`, so clearing up to the larger of this shot's and the
        // previous shot's `n` covers every node either of them could have touched.
        mark_all_unbuilt(std::max(n, last_n));
        last_n = n;

        // ---- PRE, on the manager core. Three sub-regions sharing their boundary reads, so the
        // three columns sum exactly to `ticks_pre`.
        uint64_t t0 = tick();
        union_find(h);
        uint64_t t1 = tick();
        build_tree(h);
        uint64_t t2 = tick();
        // Amendment 1 test 7, the "between the end of `tree`" half of it. Debug builds only, and a
        // register read either way.
        uint64_t writes_before_scatter = node_writes_snapshot();
        scatter(h);
        uint64_t t3 = tick();

        out.ticks_uf = t1 - t0;
        out.ticks_tree = t2 - t1;
        out.ticks_scatter = t3 - t2;
        out.ticks_pre = t3 - t0;
        check_no_manager_graph_work(writes_before_scatter);

        // The build's view of the shot. After `PRE`, because `scatter` is what sizes `refused_slot`
        // and the two `cross_index_*` arrays and a growing buffer moves; outside every timed region,
        // because in a real system these are fixed arena addresses and not work at all.
        bind_build_input(h);
        // Amendment 1 test 8: the built state does not depend on the order the leaves run in.
        check_build_independence();
        out.n_pieces = c.n_pieces;
        out.largest_piece = c.largest_piece;
        out.n_refused = (uint32_t)c.refused.size();
        out.n_fusions = c.n_fusions;
        out.tree_depth = c.depth;
        out.max_fusion_boundary = c.max_fusion_boundary;

        // ---- The job graph.
        build_jobs(out.ticks_pre);

        // ---- The schedule, driven event by event.
        bool escalated = false;
        uint32_t released = 0;
        for (size_t core = 0; core < k; core++) {
            core_free[core] = 0;
            core_last_job[core] = NONE;
        }
        manager_free = out.ticks_pre;
        manager_last_job = 0;  // PRE is job 0.
        manager_chain.clear();
        manager_chain.push_back(out.ticks_pre);

        while (!ready_heap.empty()) {
            // ---- DISPATCH, on the manager core: the ready-heap pop and the core pick.
            uint64_t d0 = tick();
            std::pop_heap(ready_heap.begin(), ready_heap.end(), [this](uint32_t a, uint32_t b) {
                const Job& x = jobs[a];
                const Job& y = jobs[b];
                if (x.priority != y.priority)
                    return x.priority < y.priority;  // Larger estimated cost first.
                return a > b;                        // Ties by job id, lower first.
            });
            uint32_t id = ready_heap.back();
            ready_heap.pop_back();
            Job& job = jobs[id];
            // `COMBINE` is the manager's own last job, not a solver job; it takes no core and is
            // charged no dispatch.
            bool on_manager = job.kind == JobKind::COMBINE;
            size_t core = 0;
            if (!on_manager) {
                if (job.core >= 0) {
                    core = (size_t)job.core;
                } else {
                    for (size_t candidate = 1; candidate < k; candidate++) {
                        if (core_free[candidate] < core_free[core])
                            core = candidate;
                    }
                }
            }
            uint64_t d1 = tick();
            uint64_t dispatch = job.charge_dispatch ? d1 - d0 : 0;
            job.dispatch_ticks = dispatch;
            out.ticks_dispatch_total += dispatch;

            // ---- The start time, and which of the three terms set it.
            uint64_t start = job.ready;
            job.predecessor = job.ready_from;
            job.determined_by = 0;
            if (!on_manager && core_free[core] > start) {
                start = core_free[core];
                job.predecessor = core_last_job[core];
                job.determined_by = 1;
            }
            if (manager_free + dispatch > start) {
                start = manager_free + dispatch;
                job.predecessor = manager_last_job;
                job.determined_by = 2;
            }
            if (job.charge_dispatch) {
                manager_chain.push_back(dispatch);
                job.manager_event = (uint32_t)manager_chain.size() - 1;
                manager_free += dispatch;
                manager_last_job = id;
            }

            // ---- The work itself.
            uint64_t w0 = tick();
            uint64_t build_ticks = 0;
            TimelineStatus status = TimelineStatus::COMPLETE;
            switch (job.kind) {
                case JobKind::LEAF: {
                    // Amendment 1 §4: a leaf is `build_piece(i)` then `solve_piece(i)`, two
                    // consecutive regions on this core, and the leaf's duration is their sum.
                    uint32_t piece = job.node;
                    build_piece(piece);
                    uint64_t b1 = tick();
                    build_ticks = b1 - w0;
                    status = fus::solve_piece(
                        shared->mwpm, decoder->horizon, c.piece_members.data() + c.piece_offsets[piece],
                        c.piece_offsets[piece + 1] - c.piece_offsets[piece]);
                    break;
                }
                case JobKind::FUSE:
                    status = fuse(h, job.node, released);
                    break;
                case JobKind::EXTRACT:
                    root_results.push_back(extract_node(job.node));
                    break;
                case JobKind::COMBINE: {
                    pm::obs_int mask = 0;
                    pm::total_weight_int weight = 0;
                    for (const pm::MatchingResult& part : root_results) {
                        mask ^= part.obs_mask;
                        weight += part.weight;
                    }
                    out.obs_mask = mask;
                    out.weight = weight;
                    break;
                }
                default:
                    break;
            }
            uint64_t w1 = tick();
            job.ticks = w1 - w0;
            job.build_ticks = build_ticks;
            job.start = start;
            job.finish = start + job.ticks;
            job.executed = true;
            if (on_manager) {
                manager_free = job.finish;
                manager_last_job = id;
            } else {
                core_free[core] = job.finish;
                core_last_job[core] = id;
            }

            switch (job.kind) {
                case JobKind::LEAF: {
                    uint64_t solve_ticks = job.ticks - job.build_ticks;
                    out.ticks_build_sum += job.build_ticks;
                    out.ticks_build_max = std::max(out.ticks_build_max, job.build_ticks);
                    out.ticks_solve_sum += solve_ticks;
                    out.ticks_solve_max = std::max(out.ticks_solve_max, solve_ticks);
                    break;
                }
                case JobKind::FUSE:
                    out.ticks_fuse_sum += job.ticks;
                    out.ticks_fuse_max = std::max(out.ticks_fuse_max, job.ticks);
                    break;
                case JobKind::EXTRACT:
                    out.ticks_extract_max = std::max(out.ticks_extract_max, job.ticks);
                    break;
                case JobKind::COMBINE:
                    out.ticks_combine = job.ticks;
                    break;
                default:
                    break;
            }

            // §11 test 1, debug builds only: the fused state stays dual-feasible.
            check_feasibility(job, status);

            if (status == TimelineStatus::TRUNCATED) {
                // Deviation 5: the shot escalates here and stops here. Its `sparse_k` is not used.
                escalated = true;
                out.escalation_site = job.kind == JobKind::LEAF ? 1 : 2;
                break;
            }
            release_dependents(id, core);
        }

        out.n_released = released;
        out.escalated = escalated ? 1 : 0;

        if (escalated) {
            // The instance keeps a shot's worth of regions and standing trees; `abandon_shot` is
            // `Mwpm::reset`, which sweeps `H`'s node vector, and is outside every region.
            abandon_shot(shared->mwpm);
            out.sparse_k = 0;
        } else {
            uint32_t combine_job = (uint32_t)jobs.size() - 1;
            out.sparse_k = jobs[combine_job].finish;
            out.cp_len = walk_critical_path(combine_job, out);
            uint64_t accounted = out.cp_manager + out.cp_build + out.cp_leaf + out.cp_fuse + out.cp_extract;
            out.cp_idle = out.sparse_k > accounted ? out.sparse_k - accounted : 0;
            // Every region was handed back by the extractions, so the cheap reset is the right one.
            reset_for_next_shot(shared->mwpm);
        }
        // No `detach` here: the mask and the dual cap are installed once at startup now, and the
        // next shot's builds write straight into the same arenas.

        // ---- FALLBACK, on its own core, on every shot so that `min(fallback, sparse_k)` exists.
        std::fill(obs_bytes.begin(), obs_bytes.end(), (uint8_t)0);
        pm::total_weight_int fallback_weight = 0;
        uint64_t f0 = tick();
        pm::decode_detection_events(*fallback_instance, shot, obs_bytes.data(), fallback_weight, false);
        uint64_t f1 = tick();
        out.fallback = f1 - f0;
        out.system = out.escalated ? out.fallback : std::min(out.fallback, out.sparse_k);
    }

    /// The observable mask of a defect's own boundary edge, taken from the ball tables the same way
    /// `BallMwpm::rebuild` does. Read outside every timed region.
    pm::obs_int boundary_obs_mask(uint64_t det) const {
        if (decoder->g_mwpm.flooder.graph.num_observables > sizeof(pm::obs_int) * 8)
            return 0;
        pm::obs_int mask = 0;
        const BallTables& tables = decoder->tables;
        for (uint64_t at = tables.bcost_mask_offsets[det]; at < tables.bcost_mask_offsets[det + 1]; at++)
            mask ^= (pm::obs_int)1 << tables.bcost_mask_ids[at];
        return mask;
    }

    /// The same, for an edge of `H`, from the ball pool entry the edge record names. Also outside
    /// every timed region: the design's edge record arrives with its mask already on it, and the
    /// build is forbidden from re-deriving it.
    pm::obs_int edge_obs_mask(uint64_t entry) const {
        if (decoder->g_mwpm.flooder.graph.num_observables > sizeof(pm::obs_int) * 8)
            return 0;
        pm::obs_int mask = 0;
        const BallTables& tables = decoder->tables;
        for (uint64_t at = tables.ball_mask_offsets[entry]; at < tables.ball_mask_offsets[entry + 1]; at++)
            mask ^= (pm::obs_int)1 << tables.ball_mask_ids[at];
        return mask;
    }

    /// How many defects the previous shot built, so that the unbuild pass covers exactly the nodes
    /// that could still be carrying adjacency.
    uint32_t last_n{0};

    /// Amendment 1 §3.1's "marks all `n_max` nodes as unbuilt". A node is built exactly when it has
    /// adjacency, so dropping the adjacency is the marking; the four vectors keep their capacity, so
    /// this frees nothing and the next shot's builds allocate nothing.
    ///
    /// Outside every timed region and charged to no core, as the amendment requires.
    void mark_all_unbuilt(uint32_t upto) {
        pm::MatchingGraph& graph = shared->mwpm.flooder.graph;
        for (uint32_t v = 0; v < upto; v++) {
            pm::DetectorNode& node = graph.nodes[v];
            node.neighbors.clear();
            node.neighbor_weights.clear();
            node.neighbor_observables.clear();
            node.neighbor_implied_weights.clear();
            built[v] = 0;
        }
    }

    /// Rebinds the build's read-only view of the shot. The pointers move when a buffer grows, so
    /// this is done once per shot rather than once per cell; it is outside every timed region.
    void bind_build_input(const BallGraph& h) {
        Cut& c = *cut;
        build_input.edges = h.edges.data();
        build_input.edge_obs = edge_obs.data();
        build_input.piece_of = c.piece_of.data();
        build_input.boundary_weight = c.bnd.data();
        build_input.boundary_is_dummy = c.via_dummy.data();
        build_input.boundary_obs = c.obs_boundary.data();
        build_input.refused_slot = refused_slot.data();
        build_input.cross_index_lo = cross_index_lo.data();
        build_input.cross_index_hi = cross_index_hi.data();
        build_input.built = built.data();
    }

    /// A root's defects are the union of its leaves' member blocks. Walked with an explicit stack
    /// rather than by recursion: the greedy minimises each fusion's boundary and not the tree's
    /// depth, so a caterpillar is a legal shape and the depth is `O(n_pieces)`.
    pm::MatchingResult extract_node(uint32_t node) {
        Cut& c = *cut;
        pm::MatchingResult result;
        walk_stack.clear();
        walk_stack.push_back(node);
        while (!walk_stack.empty()) {
            uint32_t at = walk_stack.back();
            walk_stack.pop_back();
            if (at < c.n_pieces) {
                result += fus::extract_root(
                    shared->mwpm, c.piece_members.data() + c.piece_offsets[at],
                    c.piece_offsets[at + 1] - c.piece_offsets[at]);
                continue;
            }
            walk_stack.push_back(c.tree_child1[at]);
            walk_stack.push_back(c.tree_child0[at]);
        }
        return result;
    }

    std::vector<uint32_t> walk_stack;

    /// Builds the job graph of §5: PRE, one LEAF per piece, one FUSE per internal node, one EXTRACT
    /// per root, one COMBINE. Job 0 is PRE and the last job is COMBINE.
    void build_jobs(uint64_t pre_ticks) {
        Cut& c = *cut;
        jobs.clear();
        ready_heap.clear();
        root_results.clear();

        Job pre;
        pre.kind = JobKind::PRE;
        pre.priority = UINT64_MAX;
        pre.charge_dispatch = false;
        pre.ticks = pre_ticks;
        pre.start = 0;
        pre.finish = pre_ticks;
        pre.executed = true;
        jobs.push_back(pre);

        for (uint32_t piece = 0; piece < c.n_pieces; piece++) {
            Job job;
            job.kind = JobKind::LEAF;
            job.node = piece;
            job.priority = c.cost_tab[c.piece_offsets[piece + 1] - c.piece_offsets[piece]];
            job.ready = pre_ticks;
            job.ready_from = 0;
            job.pending = 0;
            job_of_node[piece] = (uint32_t)jobs.size();
            jobs.push_back(job);
        }
        for (uint32_t node = c.n_pieces; node < c.n_tree_nodes; node++) {
            Job job;
            job.kind = JobKind::FUSE;
            job.node = node;
            job.priority = c.tree_cost[node];
            job.pending = 2;
            job_of_node[node] = (uint32_t)jobs.size();
            jobs.push_back(job);
        }
        // One extraction per root, and the combine that waits on all of them.
        uint32_t n_roots = 0;
        for (uint32_t node : c.tree_roots) {
            Job job;
            job.kind = JobKind::EXTRACT;
            job.node = node;
            // An extraction runs on its root's core the instant the root finishes, so it must win
            // any tie for that core.
            job.priority = UINT64_MAX;
            job.charge_dispatch = false;
            job.pending = 1;
            extract_of_node[node] = (uint32_t)jobs.size();
            jobs.push_back(job);
            n_roots++;
        }
        Job combine;
        combine.kind = JobKind::COMBINE;
        combine.priority = UINT64_MAX;
        combine.charge_dispatch = false;
        combine.pending = n_roots;
        combine.core = 0;
        jobs.push_back(combine);

        // A shot with no defects has no leaves and no roots: the combine is ready at once.
        if (n_roots == 0) {
            jobs.back().ready = pre_ticks;
            jobs.back().ready_from = 0;
            push_ready((uint32_t)jobs.size() - 1);
        }
        for (uint32_t piece = 0; piece < c.n_pieces; piece++)
            push_ready(job_of_node[piece]);
    }

    void push_ready(uint32_t id) {
        ready_heap.push_back(id);
        std::push_heap(ready_heap.begin(), ready_heap.end(), [this](uint32_t a, uint32_t b) {
            const Job& x = jobs[a];
            const Job& y = jobs[b];
            if (x.priority != y.priority)
                return x.priority < y.priority;
            return a > b;
        });
    }

    /// Marks `id` finished and pushes whatever that makes ready.
    void release_dependents(uint32_t id, size_t core) {
        Cut& c = *cut;
        const Job& job = jobs[id];
        uint32_t combine_job = (uint32_t)jobs.size() - 1;

        auto arrive = [&](uint32_t target, int32_t bound_core) {
            Job& next = jobs[target];
            if (job.finish > next.ready) {
                next.ready = job.finish;
                next.ready_from = id;
                // "A fusion runs on the core of its later-finishing child, which is free at that
                // instant"; an extraction runs on its root's core.
                if (bound_core >= 0)
                    next.core = bound_core;
            }
            if (--next.pending == 0)
                push_ready(target);
        };

        switch (job.kind) {
            case JobKind::LEAF:
            case JobKind::FUSE: {
                uint32_t node = job.node;
                uint32_t parent = c.tree_parent[node];
                if (parent != NONE) {
                    arrive(job_of_node[parent], (int32_t)core);
                } else {
                    arrive(extract_of_node[node], (int32_t)core);
                }
                break;
            }
            case JobKind::EXTRACT:
                arrive(combine_job, -1);
                break;
            default:
                break;
        }
    }

    /// Walks back from `COMBINE` through whatever determined each start, summing the ticks each step
    /// contributed. Returns the number of steps on the chain.
    ///
    /// Every link is tight: a dependency-determined start equals the dependency's finish, a
    /// core-determined one equals that core's free time, and a manager-determined one lands on the
    /// manager's own contiguous chain. So the chain accounts for the whole makespan and `cp_idle`
    /// comes out zero — which makes it an invariant column rather than a dead one: a non-zero
    /// `cp_idle` means the schedule and its critical path disagree.
    uint32_t walk_critical_path(uint32_t from, ShotMeasurement& out) {
        uint32_t length = 0;
        uint32_t at = from;
        while (at != NONE) {
            const Job& job = jobs[at];
            length++;
            switch (job.kind) {
                case JobKind::PRE:
                case JobKind::COMBINE:
                    out.cp_manager += job.ticks;
                    break;
                case JobKind::LEAF:
                    // Two buckets, one job: the leaf's build half and its solve half (amendment
                    // 1 §4). They sum to the leaf's contribution, so the six `cp_*` buckets still
                    // sum to `sparse_k`.
                    out.cp_build += job.build_ticks;
                    out.cp_leaf += job.ticks - job.build_ticks;
                    break;
                case JobKind::FUSE:
                    out.cp_fuse += job.ticks;
                    break;
                case JobKind::EXTRACT:
                    out.cp_extract += job.ticks;
                    break;
            }
            if (job.determined_by == 2) {
                // The manager was what held this job up, so the rest of the path is the manager's
                // own timeline: this job's dispatch, every dispatch before it, and `PRE` at the
                // bottom. `manager_chain[0]` is `PRE`'s ticks, already counted if the walk is
                // starting *at* `PRE`, which it cannot be — `PRE` has no predecessor.
                size_t top = job.manager_event != NONE ? job.manager_event : manager_chain.size() - 1;
                for (size_t i = 1; i <= top; i++) {
                    out.cp_manager += manager_chain[i];
                    length++;
                }
                out.cp_manager += manager_chain[0];
                length++;
                break;
            }
            at = job.predecessor;
            if (length > jobs.size() + manager_chain.size() + 2)
                break;  // Defensive: a cycle would be a bug, not a measurement.
        }
        return length;
    }

    /// §11 test 1. Debug builds only; compiled out entirely when `NDEBUG` is set.
    void check_feasibility(const Job& job, TimelineStatus status) {
#ifdef NDEBUG
        (void)job;
        (void)status;
#else
        if (status != TimelineStatus::COMPLETE)
            return;
        if (job.kind != JobKind::LEAF && job.kind != JobKind::FUSE)
            return;
        Cut& c = *cut;
        feasibility_scratch.clear();
        collect_defects(job.node, c, feasibility_scratch);
        const char* why = nullptr;
        if (!fus::feasibility_holds(*fusion, feasibility_scratch.data(), feasibility_scratch.size(), &why))
            feasibility_failures++;
#endif
    }

    void collect_defects(uint32_t node, Cut& c, std::vector<uint32_t>& out) {
        walk_stack.clear();
        walk_stack.push_back(node);
        while (!walk_stack.empty()) {
            uint32_t at = walk_stack.back();
            walk_stack.pop_back();
            if (at < c.n_pieces) {
                for (uint32_t i = c.piece_offsets[at]; i < c.piece_offsets[at + 1]; i++)
                    out.push_back(c.piece_members[i]);
                continue;
            }
            walk_stack.push_back(c.tree_child1[at]);
            walk_stack.push_back(c.tree_child0[at]);
        }
    }

    std::vector<uint32_t> feasibility_scratch;
    int64_t feasibility_failures{0};

    // -----------------------------------------------------------------------------------------
    // Amendment 1's tests 7, 8 and 9. All three are debug-only and compile out entirely under
    // `NDEBUG`; `fusion_lb_profiler_debug` is the target they are run on.
    // -----------------------------------------------------------------------------------------

    int64_t manager_graph_write_failures{0};
    int64_t build_independence_failures{0};
    int64_t built_before_fused_failures{0};

    uint64_t node_writes_snapshot() const {
#ifdef NDEBUG
        return 0;
#else
        return fusion->debug_node_writes;
#endif
    }

    /// Test 7. No manager-side graph work: the count of solver-node adjacency writes does not move
    /// across `scatter`, and no defect is built when `PRE` ends.
    void check_no_manager_graph_work(uint64_t before) {
#ifdef NDEBUG
        (void)before;
#else
        if (fusion->debug_node_writes != before) {
            manager_graph_write_failures++;
            return;
        }
        const pm::MatchingGraph& graph = shared->mwpm.flooder.graph;
        for (uint32_t v = 0; v < last_n; v++) {
            if (built[v] != 0 || !graph.nodes[v].neighbors.empty()) {
                manager_graph_write_failures++;
                return;
            }
        }
#endif
    }

    /// Test 8. Build independence: the leaves of one shot, run in reverse order, leave byte-identical
    /// adjacency, weights, masks, mask bytes, boundary edges and crossing-edge indices.
    ///
    /// Run before the schedule and undone afterwards, so the measured run still builds each piece
    /// exactly once, on the core and in the order the scheduler chose. It does leave the caches
    /// warmer than the scheduled build would have found them, which is a reason debug builds produce
    /// no campaign numbers and not a reason to skip the check.
    void check_build_independence() {
#ifndef NDEBUG
        Cut& c = *cut;
        if (c.n_pieces < 2)
            return;
        for (uint32_t piece = 0; piece < c.n_pieces; piece++)
            build_piece(piece);
        snapshot_built_state(forward_snapshot);
        mark_all_unbuilt(last_n);
        for (uint32_t piece = c.n_pieces; piece-- > 0;)
            build_piece(piece);
        snapshot_built_state(reverse_snapshot);
        mark_all_unbuilt(last_n);
        if (forward_snapshot != reverse_snapshot)
            build_independence_failures++;
#endif
    }

    /// Test 9. Built-before-fused: every defect of both of a fusion's children has been built by the
    /// time the fusion runs.
    void check_built_before_fused(uint32_t p) {
#ifdef NDEBUG
        (void)p;
#else
        Cut& c = *cut;
        for (uint32_t child : {c.tree_child0[p], c.tree_child1[p]}) {
            feasibility_scratch.clear();
            collect_defects(child, c, feasibility_scratch);
            for (uint32_t u : feasibility_scratch) {
                if (built[u] == 0) {
                    built_before_fused_failures++;
                    return;
                }
            }
        }
#endif
    }

#ifndef NDEBUG
    /// Every byte of the shot's built state, flattened, so that "byte-identical" is one comparison.
    std::vector<uint64_t> forward_snapshot;
    std::vector<uint64_t> reverse_snapshot;

    void snapshot_built_state(std::vector<uint64_t>& out) {
        const pm::MatchingGraph& graph = shared->mwpm.flooder.graph;
        out.clear();
        for (uint32_t v = 0; v < last_n; v++) {
            const pm::DetectorNode& node = graph.nodes[v];
            out.push_back(built[v]);
            out.push_back(node.neighbors.size());
            for (size_t k = 0; k < node.neighbors.size(); k++) {
                // The neighbour as an index rather than an address, so the snapshot says nothing
                // about where the node pool happens to live.
                out.push_back(
                    node.neighbors[k] == nullptr ? UINT64_MAX
                                                 : (uint64_t)(node.neighbors[k] - graph.nodes.data()));
                out.push_back(node.neighbor_weights[k]);
                out.push_back(node.neighbor_observables[k]);
                out.push_back(fusion->mask.get(v, (uint32_t)k));
            }
        }
        for (size_t r = 0; r < cross_index_lo.size(); r++) {
            out.push_back(cross_index_lo[r]);
            out.push_back(cross_index_hi[r]);
        }
    }
#endif

    // -----------------------------------------------------------------------------------------
    // §11 tests 3 and 4, and §8's monolithic cross-check.
    // -----------------------------------------------------------------------------------------

    /// §11 test 3: every leaf has `size <= ceil(n / k)`.
    ///
    /// §11 test 4, "cut-on-heaviest", is the claim that a cut falls on the heaviest links of a
    /// component, and it rests on three checkable facts rather than on one: the edge order really is
    /// non-decreasing in weight (amendment 1 §1 — and checking it is what makes the sorted input an
    /// assumption this binary tests rather than one it takes on trust), so the union-find really
    /// does see every lighter alternative
    /// before it refuses anything; a refused edge really does end up between two different pieces,
    /// so no refusal was wasted; and each fusion's interior list really is in non-decreasing weight,
    /// so the first crossing edge between two sides is the lightest and every other one is at least
    /// as heavy as it. Together those are §3.2's paragraph. (Checking "every crossing edge is at
    /// least as heavy as the lightest crossing edge of the same pair" on its own would be a
    /// tautology, which is why it is not what is checked.)
    void check_structure(const BallGraph& h, int64_t& threshold_violations, int64_t& heaviest_violations) {
        Cut& c = *cut;
        uint32_t n = (uint32_t)h.num_nodes();
        uint32_t threshold = (uint32_t)((n + k - 1) / k);
        if (threshold == 0)
            threshold = 1;
        for (uint32_t piece = 0; piece < c.n_pieces; piece++) {
            if (c.piece_offsets[piece + 1] - c.piece_offsets[piece] > threshold)
                threshold_violations++;
        }
        for (size_t at = 1; at < c.order.size(); at++) {
            if (h.edges[c.order[at]].w_int < h.edges[c.order[at - 1]].w_int)
                heaviest_violations++;
        }
        for (uint32_t e : c.refused) {
            if (c.piece_of[h.edges[e].i] == c.piece_of[h.edges[e].j])
                heaviest_violations++;
        }
        for (uint32_t node = c.n_pieces; node < c.n_tree_nodes; node++) {
            for (uint32_t at = c.interior_offsets[node] + 1; at < c.interior_offsets[node + 1]; at++) {
                if (h.edges[c.interior_edges[at]].w_int < h.edges[c.interior_edges[at - 1]].w_int)
                    heaviest_violations++;
            }
        }
    }

    /// §8. The monolithic solve of the whole of `H` under the same dual cap, on its own instance and
    /// outside every timed region.
    ///
    /// **The weight identity is the gate** and a mismatch throws. Observable masks are deliberately
    /// not compared: on a shot with tied optimal matchings the fused trajectory may select a
    /// different optimum, and selecting a different optimum is not a bug.
    ///
    /// The escalation decision is **counted, not asserted** — deviation 9. It is not an identity:
    /// truncation at `T` is a property of the trajectory the solver takes, not of the problem, and
    /// the fused trajectory is a different one. Both directions are safe. A shot the fused solve
    /// escalates is decoded by the fallback on `G`, which is exact. A shot it completes carries its
    /// own certificate — a feasible dual for `H` with every dual at most `T`, and a perfect matching
    /// on tight edges — and every dual being at most `T` is precisely what makes `H`'s `2T` edge set
    /// sufficient, so the answer is a global minimum-weight matching whatever the monolithic
    /// trajectory did. The two counters therefore record a difference in escalation *rate*, which is
    /// a latency question, and not a difference in correctness.
    void verify_shot(const BallGraph& h, uint64_t shot_index, const ShotMeasurement& m, CellAggregate& aggregate) {
        uint32_t n = (uint32_t)h.num_nodes();
        reference->rebuild(decoder->tables, h, decoder->arena, nullptr);
        pm::Mwpm& mono = reference->mwpm;
        mono.flooder.edge_mask = nullptr;
        mono.flooder.edge_mask_offsets = nullptr;
        mono.flooder.horizon = pm::NO_HORIZON;
        mono.flooder.dual_cap = decoder->horizon;
        mono.flooder.dual_cap_hit = false;

        all_defects.clear();
        for (uint32_t v = 0; v < n; v++)
            all_defects.push_back(v);
        TimelineStatus status =
            fus::solve_piece(mono, decoder->horizon, all_defects.data(), all_defects.size());

        bool reference_escalated = status == TimelineStatus::TRUNCATED;
        pm::total_weight_int reference_weight = 0;
        if (!reference_escalated) {
            pm::MatchingResult result = fus::extract_root(mono, all_defects.data(), all_defects.size());
            reference_weight = result.weight;
            reset_for_next_shot(mono);
        } else {
            abandon_shot(mono);
        }
        mono.flooder.dual_cap = pm::NO_HORIZON;

        if (m.escalated && !reference_escalated) {
            aggregate.escalation_divergence_fused_only++;
            if (m.escalation_site == 1)
                aggregate.escalation_divergence_fused_only_in_leaf++;
            else if (m.escalation_site == 2)
                aggregate.escalation_divergence_fused_only_in_fusion++;
            return;
        }
        if (!m.escalated && reference_escalated) {
            aggregate.escalation_divergence_monolithic_only++;
            return;
        }
        if (!reference_escalated && reference_weight != m.weight) {
            throw std::runtime_error(
                "verify: fused matching weight " + std::to_string((long long)m.weight) +
                " != monolithic " + std::to_string((long long)reference_weight) + " at shot " +
                std::to_string(shot_index));
        }
    }
};

void write_series(std::ofstream& out, const char* name, const TickSeries& series, const char* trailing) {
    std::vector<long long> samples = series.samples;
    out << "      \"" << name << "\": {\"mean\": " << fmt_g(series.mean()) << ", \"max\": " << series.max
        << ", \"p50\": " << fmt_g(percentile(samples, 0.50)) << ", \"p99\": " << fmt_g(percentile(samples, 0.99))
        << ", \"p99_9\": " << fmt_g(percentile(samples, 0.999)) << ", \"sum\": \"" << fmt_i128(series.sum)
        << "\", \"sum_sq\": \"" << fmt_i128(series.sum_sq) << "\"}" << trailing << "\n";
}

void write_histogram(std::ofstream& out, const char* name, const std::vector<uint64_t>& histogram, const char* trailing) {
    out << "      \"" << name << "\": {";
    bool first = true;
    for (size_t value = 0; value < histogram.size(); value++) {
        if (histogram[value] == 0)
            continue;
        if (!first)
            out << ", ";
        out << "\"" << value << "\": " << histogram[value];
        first = false;
    }
    out << "}" << trailing << "\n";
}

double mean_of(__int128 sum, uint64_t count) {
    return count == 0 ? 0.0 : (double)sum / (double)count;
}

void write_aggregate(std::ofstream& out, const Options& options, const std::vector<CellAggregate>& cells) {
    out << "{\n";
    out << "  \"binary\": \"fusion_lb_profiler\",\n";
    out << "  \"git_hash\": \"" << SPECMATCHING_GIT_HASH << "\",\n";
    out << "  \"stim_version\": \"" << SPECMATCHING_STIM_VERSION << "\",\n";
    out << "  \"model\": \"fusion_k_core_model\",\n";
    out << "  \"tick_unit\": \"" << hires_timer_name() << ", uncorrected; see run.log timer_overhead_ns\",\n";
    out << "  \"seed\": " << options.seed << ",\n";
    out << "  \"alpha\": " << fmt_g(options.alpha) << ",\n";
    out << "  \"verify\": " << (options.verify ? "true" : "false") << ",\n";
    out << "  \"module_means_over\": \"non-escalating shots only (deviation 5)\",\n";
    out << "  \"cells\": [\n";
    for (size_t i = 0; i < cells.size(); i++) {
        const CellAggregate& cell = cells[i];
        std::ostringstream key;
        key << "d=" << cell.distance << ",p=" << fmt_g(cell.noise) << ",T=" << fmt_g(cell.horizon)
            << ",k=" << cell.cores;
        double fallback_mean = cell.fallback.mean();
        double system_mean = cell.system.mean();
        out << "    {\n";
        out << "      \"cell\": \"" << key.str() << "\",\n";
        out << "      \"d\": " << cell.distance << ", \"p\": " << fmt_g(cell.noise)
            << ", \"T\": " << fmt_g(cell.horizon) << ", \"k\": " << cell.cores
            << ", \"alpha\": " << fmt_g(cell.alpha) << ",\n";
        out << "      \"shots\": " << cell.shots << ", \"escalations\": " << cell.escalations
            << ", \"escalations_in_leaf\": " << cell.escalations_in_leaf
            << ", \"escalations_in_fusion\": " << cell.escalations_in_fusion << ",\n";
        write_series(out, "fallback", cell.fallback, ",");
        write_series(out, "sparse_k", cell.sparse_k, ",");
        write_series(out, "system", cell.system, ",");
        uint64_t module_shots = cell.module_shots;
        out << "      \"module_means\": {\"uf\": " << fmt_g(mean_of(cell.uf_sum, module_shots))
            << ", \"tree\": " << fmt_g(mean_of(cell.tree_sum, module_shots))
            << ", \"scatter\": " << fmt_g(mean_of(cell.scatter_sum, module_shots))
            << ", \"pre\": " << fmt_g(mean_of(cell.pre_sum, module_shots))
            << ", \"dispatch_total\": " << fmt_g(mean_of(cell.dispatch_sum, module_shots))
            << ", \"build_max\": " << fmt_g(mean_of(cell.build_max_sum, module_shots))
            << ", \"build_sum\": " << fmt_g(mean_of(cell.build_sum_sum, module_shots))
            << ", \"solve_max\": " << fmt_g(mean_of(cell.solve_max_sum, module_shots))
            << ", \"solve_sum\": " << fmt_g(mean_of(cell.solve_sum_sum, module_shots))
            << ", \"fuse_max\": " << fmt_g(mean_of(cell.fuse_max_sum, module_shots))
            << ", \"fuse_sum\": " << fmt_g(mean_of(cell.fuse_sum_sum, module_shots))
            << ", \"extract_max\": " << fmt_g(mean_of(cell.extract_max_sum, module_shots))
            << ", \"combine\": " << fmt_g(mean_of(cell.combine_sum, module_shots)) << "},\n";
        out << "      \"pre_decomposition\": \"uf + tree + scatter\",\n";
        out << "      \"critical_path_means\": {\"manager\": " << fmt_g(mean_of(cell.cp_manager_sum, module_shots))
            << ", \"build\": " << fmt_g(mean_of(cell.cp_build_sum, module_shots))
            << ", \"leaf\": " << fmt_g(mean_of(cell.cp_leaf_sum, module_shots))
            << ", \"fuse\": " << fmt_g(mean_of(cell.cp_fuse_sum, module_shots))
            << ", \"extract\": " << fmt_g(mean_of(cell.cp_extract_sum, module_shots))
            << ", \"idle\": " << fmt_g(mean_of(cell.cp_idle_sum, module_shots))
            << ", \"len\": " << fmt_g(mean_of(cell.cp_len_sum, module_shots)) << "},\n";
        out << "      \"speedup_mean\": " << fmt_g(system_mean > 0 ? fallback_mean / system_mean : 0.0) << ",\n";
        write_histogram(out, "n_pieces_hist", cell.n_pieces_hist, ",");
        write_histogram(out, "tree_depth_hist", cell.tree_depth_hist, ",");
        write_histogram(out, "max_fusion_boundary_hist", cell.max_fusion_boundary_hist, ",");
        write_histogram(out, "n_released_hist", cell.n_released_hist, ",");
        out << "      \"raw_truncated_at_shot\": " << cell.raw_truncated_at_shot << ",\n";
        out << "      \"determinism_check\": " << cell.determinism_check << ",\n";
        out << "      \"verify_shots\": " << cell.verify_shots
            << ", \"escalation_divergence_fused_only\": " << cell.escalation_divergence_fused_only
            << ", \"escalation_divergence_fused_only_in_leaf\": " << cell.escalation_divergence_fused_only_in_leaf
            << ", \"escalation_divergence_fused_only_in_fusion\": " << cell.escalation_divergence_fused_only_in_fusion
            << ", \"escalation_divergence_monolithic_only\": " << cell.escalation_divergence_monolithic_only << ",\n";
        out << "      \"threshold_violations\": " << cell.threshold_violations
            << ", \"cut_on_heaviest_violations\": " << cell.cut_on_heaviest_violations
            << ", \"feasibility_violations\": " << cell.feasibility_violations
            << ", \"manager_graph_write_violations\": " << cell.manager_graph_write_violations
            << ", \"build_independence_violations\": " << cell.build_independence_violations
            << ", \"built_before_fused_violations\": " << cell.built_before_fused_violations << "\n";
        out << "    }" << (i + 1 == cells.size() ? "" : ",") << "\n";
    }
    out << "  ]\n";
    out << "}\n";
}

/// One `(d, p, T, k)` cell: warm-up, then `--shots` measured shots.
Digest run_cell(
    const Options& options,
    CellRunner& runner,
    ShotSampler& sampler,
    RowArena& arena,
    CellAggregate& aggregate,
    CellRunner::CapacityBreakdown* capacity_after_warmup,
    bool record,
    bool interactive,
    const char* label,
    size_t cell_index,
    size_t total_cells) {
    Digest digest;

    sampler.rng = std::mt19937_64(options.seed);
    size_t done = 0;
    size_t cell_total = options.warmup + options.shots;
    ShotMeasurement measurement;

    if (capacity_after_warmup != nullptr)
        *capacity_after_warmup = runner.capacity_breakdown();

    draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
    while (done < cell_total) {
        size_t batch = std::min(SHOT_BATCH, cell_total - done);
        sampler.sample(batch);
        for (size_t i = 0; i < batch; i++) {
            size_t index = done + i;
            bool measured = index >= options.warmup;
            if (measured)
                runner.run_shot<true>(sampler.shots[i], measurement);
            else
                runner.run_shot<false>(sampler.shots[i], measurement);

            if (!measured) {
                if (capacity_after_warmup != nullptr && index + 1 == options.warmup)
                    *capacity_after_warmup = runner.capacity_breakdown();
                continue;
            }
            uint64_t shot_index = (uint64_t)(index - options.warmup);

            // §11 tests 3 and 4, on every measured shot. Structural, so untimed and outside the
            // shot's regions; `arena.graph` still holds the `H` the shot decoded on.
            runner.check_structure(
                runner.decoder->arena.graph, aggregate.threshold_violations, aggregate.cut_on_heaviest_violations);
            if (options.verify) {
                runner.verify_shot(runner.decoder->arena.graph, shot_index, measurement, aggregate);
                aggregate.verify_shots++;
            }

            // §10's deterministic set: everything that is a function of `(seed, d, p, T, k, alpha)`
            // and of nothing else. Every tick column and everything derived from one is excluded,
            // because the scheduler's choices depend on measured ticks.
            digest.add(shot_index);
            digest.add(measurement.n_def);
            digest.add(measurement.n_edges);
            digest.add(measurement.n_pieces);
            digest.add(measurement.largest_piece);
            digest.add(measurement.n_refused);
            digest.add(measurement.n_fusions);
            digest.add(measurement.tree_depth);
            digest.add(measurement.max_fusion_boundary);

            if (!record)
                continue;

            Row row;
            row.shot = shot_index;
            row.m = measurement;
            arena.append(row);

            aggregate.shots++;
            aggregate.escalations += measurement.escalated;
            if (measurement.escalation_site == 1)
                aggregate.escalations_in_leaf++;
            else if (measurement.escalation_site == 2)
                aggregate.escalations_in_fusion++;
            aggregate.fallback.add(measurement.fallback);
            aggregate.system.add(measurement.system);
            if (!measurement.escalated) {
                aggregate.sparse_k.add(measurement.sparse_k);
                aggregate.module_shots++;
                aggregate.uf_sum += (__int128)measurement.ticks_uf;
                aggregate.tree_sum += (__int128)measurement.ticks_tree;
                aggregate.scatter_sum += (__int128)measurement.ticks_scatter;
                aggregate.pre_sum += (__int128)measurement.ticks_pre;
                aggregate.dispatch_sum += (__int128)measurement.ticks_dispatch_total;
                aggregate.build_max_sum += (__int128)measurement.ticks_build_max;
                aggregate.build_sum_sum += (__int128)measurement.ticks_build_sum;
                aggregate.solve_max_sum += (__int128)measurement.ticks_solve_max;
                aggregate.solve_sum_sum += (__int128)measurement.ticks_solve_sum;
                aggregate.fuse_max_sum += (__int128)measurement.ticks_fuse_max;
                aggregate.fuse_sum_sum += (__int128)measurement.ticks_fuse_sum;
                aggregate.extract_max_sum += (__int128)measurement.ticks_extract_max;
                aggregate.combine_sum += (__int128)measurement.ticks_combine;
                aggregate.cp_manager_sum += (__int128)measurement.cp_manager;
                aggregate.cp_build_sum += (__int128)measurement.cp_build;
                aggregate.cp_leaf_sum += (__int128)measurement.cp_leaf;
                aggregate.cp_fuse_sum += (__int128)measurement.cp_fuse;
                aggregate.cp_extract_sum += (__int128)measurement.cp_extract;
                aggregate.cp_idle_sum += (__int128)measurement.cp_idle;
                aggregate.cp_len_sum += (__int128)measurement.cp_len;
            }
            if (measurement.n_pieces < aggregate.n_pieces_hist.size())
                aggregate.n_pieces_hist[measurement.n_pieces]++;
            if (measurement.tree_depth < aggregate.tree_depth_hist.size())
                aggregate.tree_depth_hist[measurement.tree_depth]++;
            if (measurement.max_fusion_boundary < aggregate.max_fusion_boundary_hist.size())
                aggregate.max_fusion_boundary_hist[measurement.max_fusion_boundary]++;
            if (measurement.n_released < aggregate.n_released_hist.size())
                aggregate.n_released_hist[measurement.n_released]++;
        }
        done += batch;
        if (record && arena.rows.size() >= ROW_ARENA_FLUSH_AT)
            arena.flush();
        draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
    }
    if (record)
        arena.flush();
    if (interactive)
        std::fprintf(stderr, "\r\x1b[K");
    aggregate.feasibility_violations = runner.feasibility_failures;
    aggregate.manager_graph_write_violations = runner.manager_graph_write_failures;
    aggregate.build_independence_violations = runner.build_independence_failures;
    aggregate.built_before_fused_violations = runner.built_before_fused_failures;
    return digest;
}

/// One `(d, p)` of the grid: circuit, DEM, ball tables, the shared instance and the fallback
/// instance, then the `T` x `k` sweep over the same shots.
void run_dp_cell(
    const Options& options,
    size_t distance,
    double noise,
    RowArena& arena,
    std::vector<CellAggregate>& aggregates,
    std::ofstream& log,
    const std::string& out_dir,
    bool interactive,
    size_t& cell_index,
    size_t total_cells) {
    char setup_label[160];
    std::snprintf(setup_label, sizeof(setup_label), "d=%zu p=%g", distance, noise);
    draw_setup(interactive, setup_label, cell_index + 1, total_cells);

    // ---- §7 step 1. Circuit -> DEM, through the corpus's own generator call.
    ShotSampler sampler = ShotSampler::make(distance, distance, noise, options.seed);
    log << "\ncorpus d=" << distance << " rounds=" << distance << " p=" << noise << "\n";
    log << "  generator_call=" << generator_call(distance, distance, noise) << "\n";

    // ---- §7 step 2. Ball tables once, at `T_max = max(T list)` and `R = 2 * T_max`.
    pm::Mwpm probe = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);
    double unit = edge_weight_units(probe.flooder.graph);
    double max_horizon = *std::max_element(options.horizons.begin(), options.horizons.end());

    BallConfig ball_config;
    ball_config.T = max_horizon * unit;
    ball_config.ball.T_max = max_horizon * unit;
    ball_config.ball.R = 2.0 * ball_config.ball.T_max;
    ball_config.collect_component_stats = false;
    ball_config.collect_harvest_diagnostics = false;
    ball_config.collect_structural_counters = false;
    ball_config.verify_component_decomposition = false;
    ball_config.verify_against_g = false;
    ball_config.ball.validate();
    if (!(ball_config.ball.R >= 2.0 * ball_config.ball.T_max))
        throw std::runtime_error("the compiled ball radius is below 2 * T_max");

    BallDecoder decoder = BallDecoder::from_detector_error_model(sampler.dem, ball_config, NUM_DISTINCT_WEIGHTS);
    decoder.harvester.collect_diagnostics = false;
    decoder.harvester.use_legacy_enumeration = false;

    // The fallback core's own instance, outside the `k` and separate from `decoder.g_mwpm`, so
    // neither is warmed by the other.
    pm::Mwpm fallback_instance = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);

    const pm::MatchingGraph& g_graph = decoder.g_mwpm.flooder.graph;
    size_t n_max = g_graph.nodes.size();
    size_t num_observables = g_graph.num_observables;
    size_t k_max = options.core_counts.back();

    log << "  dem_has_negative_weights=" << (decoder.dem_has_negative_weights ? 1 : 0)
        << " edge_weight_unit=" << unit << " T_max=" << ball_config.ball.T_max << " R=" << ball_config.ball.R
        << " n_max=" << n_max << " num_observables=" << num_observables << "\n";

    // ---- §7 step 4. Everything allocated here, before the first measured shot.
    Cut cut;
    cut.configure(n_max, options.alpha);

    BallMwpm shared;
    shared.configure(num_observables, g_graph.normalising_constant);
    BallMwpm reference;
    reference.configure(num_observables, g_graph.normalising_constant);
    {
        // The shared instance is sized for the whole of `H` — the worst case is one giant component
        // — and is sized **once**, here. `rebuild` sizes the node pool from the graph it is given,
        // and after this no whole-instance rebuild happens again: per shot the graph is built piece
        // by piece inside the leaves (amendment 1 §3.1), so the `pm::Mwpm` this creates is the one
        // every shot of this `(d, p)` uses and the flooder pointers installed below stay valid.
        BallGraph full;
        full.h_to_det.resize(n_max, 0);
        shared.rebuild(decoder.tables, full, decoder.arena, nullptr);
        reference.rebuild(decoder.tables, full, decoder.arena, nullptr);
        // The sizing rebuild leaves `used_nodes = n_max`. `reference` is still rebuilt per verified
        // shot, and would otherwise clear the whole pool instead of the nodes the previous shot used.
        shared.used_nodes = 0;
        shared.grow_events = 0;
        reference.used_nodes = 0;
        reference.grow_events = 0;
    }

    fus::FusionInstance fusion;
    fusion.mask.reserve(n_max, 4 * INITIAL_EDGE_CAPACITY);

    CellRunner runner;
    runner.decoder = &decoder;
    runner.fallback_instance = &fallback_instance;
    runner.cut = &cut;
    runner.shared = &shared;
    runner.fusion = &fusion;
    runner.reference = &reference;
    runner.configure(n_max, num_observables, k_max);
    runner.refused_slot.reserve(INITIAL_EDGE_CAPACITY);
    runner.feasibility_scratch.reserve(n_max);

    log.flush();

    for (double horizon : options.horizons) {
        decoder.set_horizon(horizon * unit);
        // Attach once, here, not per shot: the mask arena and the dual cap are installed on the
        // instance that was sized above and is never replaced. `T` is the only thing that changes
        // across this loop, and it is what the cap is.
        fusion.attach(shared.mwpm, n_max, decoder.horizon);
        // `last_n` deliberately survives the `T` and `k` loops. It is what keeps "every node at or
        // above the last shot's `n` is unbuilt" true: reset it and a smaller shot would leave the
        // tail of a bigger one standing, for a later shot to build on top of.

        for (size_t k : options.core_counts) {
            cell_index++;
            char label[160];
            std::snprintf(label, sizeof(label), "d=%zu p=%g T=%g k=%zu", distance, noise, horizon, k);

            runner.k = k;
            runner.feasibility_failures = 0;
            runner.manager_graph_write_failures = 0;
            runner.build_independence_failures = 0;
            runner.built_before_fused_failures = 0;

            std::ostringstream path;
            path << out_dir << "/shots_d" << distance << "_p" << fmt_g(noise) << "_T" << fmt_g(horizon) << "_k" << k
                 << ".csv";
            std::ofstream shots_csv(path.str());
            if (!shots_csv.is_open())
                throw std::runtime_error("could not open " + path.str() + " for writing");
            shots_csv << "# d=" << distance << ",p=" << fmt_g(noise) << ",T=" << fmt_g(horizon) << ",k=" << k
                      << ",alpha=" << fmt_g(options.alpha) << ",seed=" << options.seed
                      << ",tick_unit=" << hires_timer_name() << ",uncorrected=1\n";
            shots_csv << "shot,n_def,n_edges,escalated,"
                         "n_pieces,largest_piece,n_refused,n_fusions,tree_depth,max_fusion_boundary,n_released,"
                         "ticks_uf,ticks_tree,ticks_scatter,ticks_pre,"
                         "ticks_dispatch_total,"
                         "ticks_build_max,ticks_build_sum,ticks_solve_max,ticks_solve_sum,"
                         "ticks_fuse_max,ticks_fuse_sum,ticks_extract_max,ticks_combine,"
                         "cp_manager,cp_build,cp_leaf,cp_fuse,cp_extract,cp_idle,cp_len,"
                         "ticks_sparse_k,ticks_fallback,ticks_system\n";
            arena.out = &shots_csv;
            arena.labels = CellLabels{distance, noise, horizon, k, options.alpha};
            arena.rows_this_cell = 0;
            arena.truncated_at_shot = -1;

            aggregates.emplace_back();
            CellAggregate& aggregate = aggregates.back();
            aggregate.distance = distance;
            aggregate.noise = noise;
            aggregate.horizon = horizon;
            aggregate.cores = k;
            aggregate.alpha = options.alpha;
            aggregate.reset(n_max + 2);
            aggregate.fallback.reserve(options.shots);
            aggregate.sparse_k.reserve(options.shots);
            aggregate.system.reserve(options.shots);

            log << "  cell T=" << horizon << " k=" << k << " T_weight_units=" << horizon * unit
                << " T_int=" << (int64_t)decoder.horizon << " started=" << stamp() << "\n";
            log.flush();

            CellRunner::CapacityBreakdown capacity_at_warmup_end;
            Digest digest = run_cell(
                options, runner, sampler, arena, aggregate, &capacity_at_warmup_end, true, interactive, label,
                cell_index, total_cells);
            CellRunner::CapacityBreakdown capacity_at_end = runner.capacity_breakdown();
            aggregate.raw_truncated_at_shot = arena.truncated_at_shot;

            if (options.check_determinism) {
                Digest replay = run_cell(
                    options, runner, sampler, arena, aggregate, nullptr, false, interactive, label, cell_index,
                    total_cells);
                aggregate.determinism_check = replay.value == digest.value ? 1 : 0;
            }

            log << "  cell T=" << horizon << " k=" << k << " finished=" << stamp() << " rows=" << arena.rows_this_cell
                << " escalations=" << aggregate.escalations
                << " no_growth_after_warmup=" << (capacity_at_end == capacity_at_warmup_end ? 1 : 0)
                << " (grown: input=" << (capacity_at_end.input - capacity_at_warmup_end.input)
                << " solver=" << (capacity_at_end.solver - capacity_at_warmup_end.solver)
                << " own=" << (capacity_at_end.own - capacity_at_warmup_end.own) << ")"
                << " threshold_violations=" << aggregate.threshold_violations
                << " cut_on_heaviest_violations=" << aggregate.cut_on_heaviest_violations
                << " feasibility_violations=" << aggregate.feasibility_violations
                << " manager_graph_write_violations=" << aggregate.manager_graph_write_violations
                << " build_independence_violations=" << aggregate.build_independence_violations
                << " built_before_fused_violations=" << aggregate.built_before_fused_violations;
            if (options.verify) {
                log << " verify_shots=" << aggregate.verify_shots
                    << " escalation_divergence_fused_only=" << aggregate.escalation_divergence_fused_only
                    << " (leaf=" << aggregate.escalation_divergence_fused_only_in_leaf
                    << ",fusion=" << aggregate.escalation_divergence_fused_only_in_fusion << ")"
                    << " escalation_divergence_monolithic_only=" << aggregate.escalation_divergence_monolithic_only;
            }
            if (options.check_determinism)
                log << " determinism_check=" << (aggregate.determinism_check == 1 ? "pass" : "FAIL");
            if (arena.truncated_at_shot >= 0)
                log << " raw_truncated_at_shot=" << arena.truncated_at_shot;
            log << "\n";
            log.flush();

            arena.out = nullptr;
            shots_csv.close();

            if (options.check_determinism && aggregate.determinism_check != 1)
                throw std::runtime_error(
                    "--check-determinism: the structural columns differ between two runs of " + std::string(label));
            if (aggregate.threshold_violations != 0)
                throw std::runtime_error("test 3 (threshold): a piece exceeded ceil(n/k) at " + std::string(label));
            if (aggregate.cut_on_heaviest_violations != 0)
                throw std::runtime_error("test 4 (cut-on-heaviest) failed at " + std::string(label));
            if (aggregate.feasibility_violations != 0)
                throw std::runtime_error("test 1 (feasibility invariant) failed at " + std::string(label));
            if (aggregate.manager_graph_write_violations != 0)
                throw std::runtime_error(
                    "amendment test 7 (no manager-side graph work) failed at " + std::string(label));
            if (aggregate.build_independence_violations != 0)
                throw std::runtime_error("amendment test 8 (build independence) failed at " + std::string(label));
            if (aggregate.built_before_fused_violations != 0)
                throw std::runtime_error("amendment test 9 (built-before-fused) failed at " + std::string(label));
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }

    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }

    std::ofstream log(options.out_dir + "/run.log");
    if (!log.is_open()) {
        std::cerr << "error: could not open " << options.out_dir << "/run.log\n";
        return 1;
    }

    RowArena arena;
    arena.configure();

    // The startup calibration, before any decoding, so it measures the clock rather than the clock
    // plus whatever the decoder left in the caches.
    double overhead_ns = hires_timer_overhead_ns();
    bool thread_scoped = hires_timer_is_thread_scoped();

    log << "fusion_lb_profiler - divide-and-conquer sparse blossom on H with a k-core list scheduler\n";
    log << "git_hash=" << SPECMATCHING_GIT_HASH << "\n";
    log << "build_flags=" << SPECMATCHING_BUILD_FLAGS << "\n";
    log << "stim_version=" << SPECMATCHING_STIM_VERSION << "\n";
#if defined(__APPLE__)
    log << "platform=macos\n";
#else
    log << "platform=linux\n";
#endif
    log << "hires_timer_name=" << hires_timer_name() << "\n";
    log << "timer_is_thread_scoped=" << (thread_scoped ? 1 : 0) << "\n";
    log << "timer_overhead_ns=" << overhead_ns << " (10^6 back-to-back hires_now_ns() calls)\n";
    log << "rows_are_uncorrected=1 (every ticks_* column is a raw delta of this clock; the analysis"
           " subtracts timer_overhead_ns)\n";
    log << "threads=1 (no thread is ever spawned)\n";
    log << "fusion_k_core_model=1: one shared pm::Mwpm holds the shot's whole H; every unit of work"
           " runs on the profiling thread and is timed, and a simulated list scheduler places those"
           " measurements on k solver cores plus one manager core plus one fallback core. sparse_k is"
           " the makespan of that schedule.\n";
    log << "shared_instance=1, transfer_cost_modelled=0: nothing is copied between cores and the"
           " cache-miss cost of touching another core's data is NOT modelled.\n";
    log << "warm_cache_caveat=1: a unit of work executed immediately after another sees a warmer"
           " cache than a real core would; this travels with every number in every CSV and in"
           " agg.json.\n";
    log << "timed_region_input=compute_seeded_detection_events + build_ball_graph are the INPUT and"
           " are outside every timed region; so are the edge-list sort, the per-edge and per-defect"
           " observable lookups, marking the shared instance's nodes unbuilt, the instance reset"
           " (reset_for_next_shot / abandon_shot) and the --verify cross-check\n";
    log << "edge_order=sorted_by_hardware (std::sort in simulator, untimed): the edge list reaches"
           " the solver in non-decreasing weight, ties by (u, v), because the sort happens in the"
           " hardware that produces it. No bucket arena exists and there is no ticks_bucket column.\n";
    log << "graph_build=per_piece_in_leaf (no whole-instance rebuild per shot): each leaf writes the"
           " adjacency, weights, observable masks, mask bytes and boundary edge of its own defects,"
           " on its own core, as the first of its two timed regions. The manager used to do all of"
           " it, serially, inside scatter.\n";
    log << "adjacency_order=insertion_sort_per_node: a piece's edge list arrives in weight order and"
           " adjacency must be written in ascending neighbour index, so each defect's entries are"
           " gathered into the edge mask's own offsets and insertion-sorted there. Degrees in H are"
           " small and the design offers the choice; the gather slices are per node, so no two"
           " pieces address the same word.\n";
    log << "node_pool_clearing=never_charged_to_any_core: a node is built exactly when it has"
           " adjacency, so dropping the adjacency is what marks it unbuilt. That pass sits beside"
           " reset_for_next_shot in the untimed run of the shot loop - at the top of the next shot"
           " rather than the bottom of this one, which is the same place in the loop and makes the"
           " first shot's state explicit. The four vectors keep their capacity, so it frees nothing"
           " and the next shot's builds allocate nothing.\n";
    log << "attach=once_at_startup: the mask arena and the dual cap are installed on the shared"
           " instance once per (d, p, T), after the one sizing rebuild at n_max, and not per shot."
           " Nothing replaces the pm::Mwpm any more, so nothing invalidates the flooder's pointers;"
           " the mask offsets are re-derived per shot from the degrees the builds are about to"
           " write, and the arena pointers are re-installed with them in case a buffer grew.\n";
    log << "timed_region_pre=uf + tree + scatter, three sub-regions sharing their boundary reads so"
           " the three columns sum exactly to ticks_pre\n";
    log << "timed_region_uf=bounded union-find at S = ceil(n/k) with the per-defect boundary"
           " weights, then the pieces, the per-piece defect blocks and the quotient graph\n";
    log << "timed_region_tree=cheapest-boundary-first fusion tree, then each fusion's interior"
           " crossing-edge list\n";
    log << "timed_region_scatter=the sorted edge list grouped by piece (counting sort; an interior"
           " edge once, a crossing edge once under each endpoint's piece), the per-defect adjacency"
           " degrees, the per-defect crossing-edge lists, refused_slot, the mask/gather offsets and"
           " the tuf reset. O(n + m) integer work in fixed arenas: it reads no ball table and"
           " touches no solver node, which amendment test 7 checks outright.\n";
    log << "timed_region_dispatch=the ready-heap pop and the core pick, charged to the manager core,"
           " once per LEAF and per FUSE. EXTRACT and COMBINE are charged no dispatch: an extraction"
           " runs on its root's core the instant the root finishes, and the combine is the manager's"
           " own last job.\n";
    log << "timed_region_build=build_piece(i): write the adjacency, weights, observable masks, mask"
           " bytes and boundary edge of this piece's defects and of no others. The far node of a"
           " crossing edge is an address into the shared node array and is valid before that node is"
           " built, and the crossing edge is masked, so the flooder never follows it before the"
           " fusion that unmasks it - by which time both sides are built.\n";
    log << "timed_region_leaf=solve_piece(i): inject the piece's defects at the shared clock, then"
           " run the flooder until every region of the piece is matched or a dual cap fires. The"
           " build above and this are two consecutive regions on one core and the leaf's duration is"
           " their sum; ticks_build_* and ticks_solve_* report them separately, as do cp_build and"
           " cp_leaf.\n";
    log << "timed_region_fuse=fuse(p): unmask the crossing edges that become interior, raise the"
           " dummies they stood for, release the regions those dummies held matched, run the"
           " flooder\n";
    log << "timed_region_extract=extract(r): shatter and accumulate the observable mask and weight"
           " of the root's matched pairs\n";
    log << "timed_region_combine=XOR of the root observable masks, sum of the root weights\n";
    log << "timed_region_fallback=pm::decode_detection_events on G from the FULL syndrome, on its"
           " own instance, run on every shot (not only escalating ones)\n";
    log << "sparse_k=finish time of COMBINE in the simulated schedule\n";
    log << "system=escalated ? fallback : min(fallback, sparse_k); nothing else in this artifact is"
           " called system\n";
    log << "horizon_definition=a cap on DUALS, not on the clock: the run escalates when a growing"
           " region would take some member defect's dual past T. Enforced per region, because pieces"
           " share one clock and released regions resume with radii they earned earlier, so the"
           " clock tracks no single dual. inner_max(R) is the largest frozen (blossom-nested) part"
           " among R's member defects and the cap is radius(R) = T - inner_max(R).\n";
    log << "decoder_change_1=GraphFlooder::edge_mask/edge_mask_offsets - a byte per directed"
           " adjacency entry; a masked half-edge generates no collision event. nullptr on every"
           " instance but this binary's shared one, and the neighbour scan is templated on its"
           " presence, so the stock decode path is the code it was.\n";
    log << "decoder_change_2=GraphFlooder::dual_cap/dual_cap_hit and"
           " GraphFlooder::schedule_dual_cap_event - the per-region horizon above. The cap event"
           " rides the region's shrink_event_tracker, which is idle for exactly as long as the"
           " region is growing, so GraphFillRegion is not widened. NO_HORIZON on every instance but"
           " this binary's.\n";
    log << "decoder_otherwise_unchanged=1 (nothing else under src/specmatching/ is edited)\n";
    log << "deviation_1=file lives in benchmarks/spec_matching/, not benchmarks/two_phase/; the"
           " generator call is this corpus's own (rotated_memory_x, three noise parameters); --T is"
           " in multiples of one lattice edge weight, not raw DEM float units\n";
    log << "deviation_2=the design's bucketing pass (§3.1) is gone: amendment 1 has the edge list"
           " arrive already sorted, so the simulator's stand-in std::sort is untimed and w_1 and the"
           " bucket count are no longer used anywhere. The refusal property needs only non-decreasing"
           " weight, which the sorted input gives, and (w, u, v) is a strict total order on H's edges"
           " so the permutation is the one bucketing would have produced.\n";
    log << "deviation_3=the quotient graph is keyed on the pieces' FINAL ids, built from the refused"
           " list at the end of the same timed uf region, because a union-find root at the moment of"
           " refusal is not final\n";
    log << "deviation_4=a fusion's interior crossing-edge list is assigned by lowest common ancestor"
           " once the tree is built, rather than accumulated during the greedy's adjacency merge\n";
    log << "deviation_5=an escalating shot stops at the job that truncated; its remaining jobs are"
           " NOT executed, because a truncated node leaves a primal-dual state that fusion is not"
           " defined on. agg.json takes every ticks_* and cp_* mean over NON-ESCALATING shots only,"
           " and sparse_k is 0 on an escalating row.\n";
    log << "deviation_6=per-node adjacency inside a build is ordered by insertion sort; see"
           " adjacency_order above\n";
    log << "deviation_7=a defect left with no boundary at all has its boundary edge raised to 2T + 2"
           " (unreachable, since no dual exceeds T) rather than erased\n";
    log << "deviation_8=the dummy boundary weight is floor(w/2) rounded DOWN TO EVEN. Two regions"
           " growing towards each other meet at (w - y1 - y2) >> 1, which is exact only when the"
           " summed duals are even; a monolithic run gets that for free, two independently solved"
           " pieces do not, and without the rounding --verify fails by exactly 1 unit of weight on"
           " some shots. A released region was matched to a dummy and a boundary match is tight, so"
           " its dual IS the dummy weight; making every dummy even makes every released dual even."
           " Rounding down only strengthens bnd(u) + bnd(v) <= w. The other half of the same rule:"
           " growth must also BEGIN on an even clock, so the shared clock is nudged to the next even"
           " tick before a piece is injected and before a fusion releases anything. Both points sit"
           " between solves, where the queue is drained, so the nudge disturbs no scheduled event.\n";
    log << "balancer=cut inside the union-find: an edge whose union would take a piece past"
           " S = ceil(n/k) is refused and becomes a crossing edge. Edges are visited in"
           " non-decreasing weight, so a cut falls on the heaviest links of a component.\n";
    log << "fusion_tree=cheapest-boundary-first greedy (Huffman with the crossing count as the merge"
           " cost), min-heap keyed on (crossing_count, min_id, max_id), stale entries discarded on"
           " pop. It minimises the boundary of each individual fusion, not the depth of the tree;"
           " tree_depth is logged so the trade can be seen.\n";
    log << "scheduler=list scheduler over k solver cores. Ready jobs are prioritised by estimated"
           " cost, larger first (cost[size] for a leaf, the crossing count for a fusion), ties by job"
           " id. A fusion runs on the core of its later-finishing child; an extraction on its root's"
           " core. A job starts at max(ready time, core free time, manager free time + dispatch"
           " ticks). Execution order on the thread is the dispatch order, which is a topological"
           " order of the fusion tree.\n";
    log << "amendment_tests=7 (no manager-side graph work: the count of solver-node adjacency writes"
           " does not move across scatter, and no defect is built when PRE ends), 8 (build"
           " independence: the leaves run in reverse order leave byte-identical adjacency, weights,"
           " masks, mask bytes and crossing-edge indices), 9 (built-before-fused: every defect of"
           " both children of a fusion is built when the fusion runs). All three are debug-only and"
           " are run by the fusion_lb_profiler_debug target; a non-zero count aborts the cell.\n";
    log << "critical_path=walked back from COMBINE through whatever determined each start: the"
           " dependency, the core, or the manager. A manager-determined start chains into the"
           " manager's OWN contiguous timeline (PRE, then every dispatch up to this job's), not into"
           " the job that was dispatched last - that job is still running on a solver core. Every"
           " link is therefore tight and cp_idle is 0; a non-zero cp_idle means the schedule and its"
           " critical path disagree, which is a bug rather than a measurement.\n";
    log << "alpha=" << options.alpha
        << " (leaf cost exponent for the scheduler's priority; at alpha = 1 the leaf cost is the"
           " defect count and no floating-point operation is on the critical path)\n";
    log << "verify=" << (options.verify ? 1 : 0)
        << " (per shot, the monolithic solve of the whole of H under the same dual cap, on its own"
           " instance, untimed. The total matching WEIGHT is the gate and a mismatch aborts the run."
           " Observable masks are NOT compared: on a shot with tied optimal matchings the fused"
           " trajectory may select a different optimum.)\n";
    log << "deviation_9=the escalation decision is COUNTED, not asserted. The design asks --verify to"
           " check that the fused and monolithic solves both escalate or neither does; that is not a"
           " theorem and it fails on real shots, because truncation at T is a property of the"
           " TRAJECTORY and a fused solve takes a different one (its pieces reach a dual optimum of"
           " their own, and a released region resumes from a radius the monolithic run never gave"
           " it). Both directions are safe: an escalated shot goes to the exact fallback on G, and a"
           " completed one carries a feasible dual for H with every dual <= T plus a perfect matching"
           " on tight edges, which is precisely the certificate that makes H's 2T edge set sufficient"
           " - so it is a global MWPM whatever the monolithic trajectory did. agg.json reports"
           " escalation_divergence_fused_only (split leaf/fusion) and"
           " escalation_divergence_monolithic_only per cell; they are an escalation-RATE difference,"
           " which is a campaign question.\n";
    log << "check_determinism=" << (options.check_determinism ? 1 : 0)
        << " (diffs the structural columns only: n_def, n_edges, n_pieces, largest_piece, n_refused,"
           " n_fusions, tree_depth, max_fusion_boundary. Every tick column and everything derived"
           " from one is excluded, because the scheduler's choices depend on measured ticks.)\n";
    log << "seed=" << options.seed << " (shot i of (d, p) is identical across the T and k sweeps)\n";
    log << "shots_per_cell=" << options.shots << " (after warm-up)\n";
    log << "warmup_shots_per_cell=" << options.warmup << " (run and discarded)\n";
    log << "warmup_is_not_profiled=1: a warm-up shot runs the identical workload and reads no clock"
           " at all; run_shot is templated on the switch.\n";
    log << "T_unit=one lattice edge weight (the median discretised edge of G, in DEM float units);"
           " each cell below records T, T_weight_units and T_int\n";

    if (!thread_scoped) {
        log << "WARNING: no thread-scoped clock was available, so every tick on this run includes"
               " time the thread spent off the CPU and is an upper bound only.\n";
        std::fprintf(
            stderr,
            "\n  WARNING: latencies are being measured with a WALL CLOCK, backend '%s'.\n"
            "  Every tick includes time the thread spent off the CPU and is an upper bound only.\n\n",
            hires_timer_name());
        std::fflush(stderr);
    }
    log.flush();

    bool interactive = stderr_is_terminal();
    size_t total_cells = options.distances.size() * options.error_rates.size() * options.horizons.size() *
                         options.core_counts.size();
    size_t cell_index = 0;
    bool ok = true;
    std::vector<CellAggregate> aggregates;
    aggregates.reserve(total_cells);

    for (size_t distance : options.distances) {
        if (!ok)
            break;
        for (double noise : options.error_rates) {
            try {
                run_dp_cell(
                    options, distance, noise, arena, aggregates, log, options.out_dir, interactive, cell_index,
                    total_cells);
            } catch (const std::exception& error) {
                std::fprintf(stderr, "\n");
                std::cerr << "error at d=" << distance << " p=" << noise << ": " << error.what() << "\n";
                log << "  ERROR=" << error.what() << "\n";
                ok = false;
                break;
            }
        }
    }

    {
        std::ofstream agg(options.out_dir + "/agg.json");
        if (agg.is_open()) {
            write_aggregate(agg, options, aggregates);
        } else {
            std::cerr << "warning: could not open " << options.out_dir << "/agg.json for writing\n";
            log << "WARNING: could not open agg.json for writing\n";
        }
    }

    log << "\nfinished=" << stamp() << "\n";
    log << "rows_total=" << arena.rows_total << "\n";
    log << (ok ? "status=ok\n" : "status=incomplete (see ERROR above)\n");

    std::printf("wrote %s/{shots_*.csv,agg.json,run.log}\n", options.out_dir.c_str());
    return ok ? 0 : 1;
}
