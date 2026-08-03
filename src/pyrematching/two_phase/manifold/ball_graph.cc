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
    BallGraphTiming* timing) {
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
    for (uint32_t i = 0; i < graph.h_to_det.size(); i++) {
        uint64_t det = graph.h_to_det[i];

        if (mode == BallGraphBuildMode::SCAN) {
            uint64_t end = tables.ball_end(det);
            for (uint64_t e = tables.ball_begin(det); e < end; e++) {
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
                    if ((pm::cumulative_time_int)w_int > two_t)
                        continue;
                    graph.edges.push_back(BallGraphEdge{i, arena.det_to_h[target], w_int, entry});
                }
            }
        }

        if (tables.has_bcost[det] &&
            (pm::cumulative_time_int)tables.bcost_w_int[det] <= (pm::cumulative_time_int)horizon)
            graph.boundary_edges.push_back(BallBoundaryEdge{i, tables.bcost_w_int[det], det});
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
