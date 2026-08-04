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

#include "pyrematching/two_phase/manifold/ball_graph.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>

#include "pyrematching/two_phase/perf/two_phase_profile.h"

namespace pm {
namespace two_phase {

namespace {

constexpr uint32_t NOT_A_DEFECT = std::numeric_limits<uint32_t>::max();

/// The sum of every buffer's capacity. Comparing it before and after a build is a cheap, allocator
/// -independent way to notice that the arena grew.
uint64_t total_capacity(const BallGraphArena& arena) {
    return arena.graph.h_to_det.capacity() + arena.graph.edges.capacity() + arena.graph.boundary_edges.capacity() +
           arena.det_to_h.capacity() + arena.syndrome_words.capacity() + arena.touched_words.capacity();
}

/// Bytes of the observable id list behind one ball entry, plus the two CSR offsets that bound it.
/// Physically these are read by `BallMwpm::rebuild`, not by the loop below; they are charged to the
/// intersection because the pair is what selects them, and because the hardware stage the counter
/// is sizing does the selection and the fetch together.
uint64_t mask_bytes(const BallTables& tables, uint64_t entry) {
    using namespace ball_element_bytes;
    return 2 * MASK_OFFSET + (tables.ball_mask_offsets[entry + 1] - tables.ball_mask_offsets[entry]) * MASK_ID;
}

uint64_t boundary_mask_bytes(const BallTables& tables, uint64_t det) {
    using namespace ball_element_bytes;
    return 2 * BCOST_MASK_OFFSET +
           (tables.bcost_mask_offsets[det + 1] - tables.bcost_mask_offsets[det]) * BCOST_MASK_ID;
}

/// What `BITSET` would have read for `det`, without running it: the node's bitset window, clamped
/// exactly where the loop clamps it — a window whose tail runs past the syndrome bitset is cut
/// short there.
uint64_t bitset_words_counterfactual(const BallTables& tables, uint64_t det, size_t syndrome_words) {
    uint64_t base = tables.ball_word_base[det];
    if (base >= syndrome_words)
        return 0;
    return std::min(tables.word_len(det), (uint64_t)syndrome_words - base);
}

/// What `SCAN` would have walked for `det` at this horizon, without running it. The entries are
/// sorted by `(w_int, target)`, so the walk stops at the first one past `2T` — and that entry's
/// weight is still read, which is the `+ 1`.
uint64_t scan_entries_counterfactual(const BallTables& tables, uint64_t det, pm::cumulative_time_int two_t) {
    auto first = tables.ball_w_int.begin() + (ptrdiff_t)tables.ball_begin(det);
    auto last = tables.ball_w_int.begin() + (ptrdiff_t)tables.ball_end(det);
    auto stop = std::upper_bound(first, last, two_t, [](pm::cumulative_time_int bound, pm::weight_int w) {
        return bound < (pm::cumulative_time_int)w;
    });
    uint64_t walked = (uint64_t)(stop - first);
    return stop == last ? walked : walked + 1;
}

}  // namespace

void BallGraphArena::reset_for_graph(size_t num_detector_nodes) {
    graph.clear();
    det_to_h.assign(num_detector_nodes, NOT_A_DEFECT);
    syndrome_words.assign((num_detector_nodes + 63) / 64, 0);
    touched_words.clear();
    grow_events = 0;
}

void build_ball_graph(
    const BallTables& tables,
    const std::vector<uint64_t>& seeded_dets,
    horizon_int horizon,
    BallGraphArena& arena,
    BallGraphBuildMode mode,
    BallGraphTiming* timing,
    BallGraphCounts* counts) {
    if (horizon != pm::NO_HORIZON && horizon > tables.t_max_int)
        throw std::invalid_argument(
            "The requested horizon exceeds the compiled BallParams::T_max; recompile the ball tables. Decoding "
            "beyond T_max would silently drop reachable pairs.");
    if (arena.det_to_h.size() != tables.num_nodes)
        arena.reset_for_graph(tables.num_nodes);

    uint64_t capacity_before = total_capacity(arena);
    BallGraph& graph = arena.graph;
    graph.clear();

    // The horizon is finite in every truncated run; `NO_HORIZON` would ask for the complete defect
    // metric, which the ball tables do not hold and which `R >= 2 * T_max` was never sized for.
    if (horizon == pm::NO_HORIZON)
        throw std::invalid_argument(
            "The ball graph front end needs a finite horizon: at T = infinity `H` is the complete defect graph, "
            "which no ball radius covers. Decode that case on G with M1's path.");
    pm::cumulative_time_int two_t = 2 * (pm::cumulative_time_int)horizon;

    HiResTimer timer;
    if (timing != nullptr)
        timer.start();

    // ---- Nodes, and the two lookup structures the intersection reads.
    graph.h_to_det.reserve(seeded_dets.size());
    for (uint64_t det : seeded_dets) {
        if (det >= tables.num_nodes)
            throw std::invalid_argument("Detection event index is outside the detector graph.");
        assert(
            (graph.h_to_det.empty() || det > graph.h_to_det.back()) &&
            "seeded detection events must be sorted ascending and deduplicated");
        arena.det_to_h[det] = (uint32_t)graph.h_to_det.size();
        graph.h_to_det.push_back(det);

        size_t word = det / 64;
        if (arena.syndrome_words[word] == 0)
            arena.touched_words.push_back((uint32_t)word);
        arena.syndrome_words[word] |= (uint64_t)1 << (det % 64);
    }

    // ---- Edges. Each pair is emitted once, from its lower-id endpoint.
    //
    // The `structural_*` locals are the §M2 hardware-budget counters. Everything that is a plain
    // register increment is accumulated unconditionally and stored out once at the end, so the
    // no-profile path pays a loop-carried add and no memory traffic; anything that needs a table
    // read the decode itself does not do sits behind `counts != nullptr`.
    uint64_t structural_scan_entries = 0;
    uint64_t structural_scan_words = 0;
    uint64_t structural_hit_bytes = 0;
    uint64_t structural_other_mode_bytes = 0;
    uint64_t structural_edges = 0;
    uint64_t structural_boundary_edges = 0;
    for (uint32_t i = 0; i < graph.h_to_det.size(); i++) {
        uint64_t det = graph.h_to_det[i];

        if (mode == BallGraphBuildMode::SCAN) {
            uint64_t end = tables.ball_end(det);
            for (uint64_t e = tables.ball_begin(det); e < end; e++) {
                structural_scan_entries++;
                // Entries are sorted by `(w_int, target)`, so the first entry past `2T` ends the
                // scan: everything after it is unreachable by §M2.0, not merely unlikely.
                if ((pm::cumulative_time_int)tables.ball_w_int[e] > two_t)
                    break;
                uint32_t target = tables.ball_target[e];
                if (target <= det)
                    continue;
                uint32_t j = arena.det_to_h[target];
                if (j == NOT_A_DEFECT)
                    continue;
                graph.edges.push_back(BallGraphEdge{i, j, tables.ball_w_int[e], e});
                structural_edges++;
                if (counts != nullptr)
                    structural_hit_bytes += mask_bytes(tables, e);
            }
        } else {
            uint64_t words_begin = tables.ball_word_offsets[det];
            uint64_t words_end = tables.ball_word_offsets[det + 1];
            uint32_t base = tables.ball_word_base[det];
            uint64_t entry_begin = tables.ball_begin(det);
            for (uint64_t k = words_begin; k < words_end; k++) {
                size_t syndrome_word = base + (k - words_begin);
                if (syndrome_word >= arena.syndrome_words.size())
                    break;
                structural_scan_words++;
                uint64_t hits = tables.ball_words[k] & arena.syndrome_words[syndrome_word];
                while (hits != 0) {
                    uint64_t bit = hits & (~hits + 1);
                    int index = __builtin_ctzll(hits);
                    hits ^= bit;
                    uint32_t target = (uint32_t)(syndrome_word * 64 + index);
                    if (target <= det)
                        continue;
                    // Rank within the node's bitset is the index into the ascending-target entry
                    // pool, which is how the bitset path names the same canonical entry SCAN does.
                    uint32_t rank =
                        tables.ball_word_rank[k] + (uint32_t)__builtin_popcountll(tables.ball_words[k] & (bit - 1));
                    uint64_t entry = tables.ball_entry_by_rank[entry_begin + rank];
                    pm::weight_int w_int = tables.ball_w_int[entry];
                    // The rank indirection and the weight are read for every candidate, including
                    // the ones the `2T` test then rejects: the bitset carries neither, so there is
                    // no way to apply the test without fetching them first.
                    if (counts != nullptr)
                        structural_hit_bytes += ball_element_bytes::WORD_RANK + ball_element_bytes::ENTRY_BY_RANK +
                                                ball_element_bytes::WEIGHT;
                    if ((pm::cumulative_time_int)w_int > two_t)
                        continue;
                    graph.edges.push_back(BallGraphEdge{i, arena.det_to_h[target], w_int, entry});
                    structural_edges++;
                    if (counts != nullptr)
                        structural_hit_bytes += mask_bytes(tables, entry);
                }
            }
        }

        bool has_boundary = tables.has_bcost[det] != 0;
        bool boundary_within_horizon =
            has_boundary && (pm::cumulative_time_int)tables.bcost_w_int[det] <= (pm::cumulative_time_int)horizon;
        if (boundary_within_horizon) {
            graph.boundary_edges.push_back(BallBoundaryEdge{i, tables.bcost_w_int[det], det});
            structural_boundary_edges++;
        }
        if (counts != nullptr) {
            structural_hit_bytes += ball_element_bytes::HAS_BCOST;
            if (has_boundary)
                structural_hit_bytes += ball_element_bytes::BCOST_WEIGHT;
            if (boundary_within_horizon)
                structural_hit_bytes += boundary_mask_bytes(tables, det);
            structural_other_mode_bytes +=
                mode == BallGraphBuildMode::SCAN
                    ? bitset_words_counterfactual(tables, det, arena.syndrome_words.size()) * ball_element_bytes::WORD
                    : scan_entries_counterfactual(tables, det, two_t) *
                          (ball_element_bytes::TARGET + ball_element_bytes::WEIGHT);
        }
    }

    // The edge records emitted are exactly what `H` ends up holding: `graph.edges` is cleared at
    // the top and only ever appended to, and each undirected pair is appended once, from its
    // lower-id endpoint. That is what lets the aggregate read `hbld_edges_written` off `H`'s own
    // sizes and still call it a count of writes.
    assert(structural_edges == graph.edges.size() && "an undirected pair was emitted more than once");
    assert(structural_boundary_edges == graph.boundary_edges.size());

    if (counts != nullptr) {
        counts->isect_scan_bytes = structural_scan_entries * (ball_element_bytes::TARGET + ball_element_bytes::WEIGHT) +
                                   structural_scan_words * ball_element_bytes::WORD;
        counts->isect_hit_bytes = structural_hit_bytes;
        counts->isect_scan_bytes_other_mode = structural_other_mode_bytes;
        counts->edges_written = structural_edges;
        counts->boundary_edges_written = structural_boundary_edges;
    }

    if (timing != nullptr) {
        timing->intersect_ns = timer.elapsed_ns();
        timer.start();
    }

    // SCAN emits a node's hits in weight order and BITSET in target order. Sorting here is what
    // makes the two modes produce the *same* `H`, and what keeps `H` a function of the shot alone.
    std::sort(graph.edges.begin(), graph.edges.end(), [](const BallGraphEdge& a, const BallGraphEdge& b) {
        return a.i != b.i ? a.i < b.i : a.j < b.j;
    });

    // ---- Undo the syndrome bitset through the touched list; `det_to_h` through the node list.
    for (uint32_t word : arena.touched_words)
        arena.syndrome_words[word] = 0;
    arena.touched_words.clear();
    for (uint64_t det : graph.h_to_det)
        arena.det_to_h[det] = NOT_A_DEFECT;

    if (total_capacity(arena) != capacity_before)
        arena.grow_events++;
    if (timing != nullptr)
        timing->finalize_ns = timer.elapsed_ns();
}

}  // namespace two_phase
}  // namespace pm
