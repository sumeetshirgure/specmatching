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

#ifndef PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_TABLES_H
#define PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_TABLES_H

#include <algorithm>
#include <cstdint>
#include <vector>

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "pyrematching/two_phase/manifold/ball_params.h"
#include "pyrematching/two_phase/truncation/horizon.h"

namespace pm {
namespace two_phase {

/// What compiling the tables cost and what it produced. Reported in the M2 exit artifact: the table
/// size is the headline number for the local-memory hardware argument.
struct BallStats {
    uint64_t num_nodes{0};
    uint64_t total_entries{0};
    double mean_ball_size{0};
    uint64_t max_ball_size{0};
    /// Entry counts per shell, summed over nodes. `shell_entry_counts.size() == num_shells`.
    std::vector<uint64_t> shell_entry_counts;

    uint64_t nodes_with_boundary{0};

    double mean_ball_word_len{0};
    uint64_t max_ball_word_len{0};

    /// Ball entries whose observable set changes when the Dijkstra tie-break is perturbed, i.e.
    /// pairs with two shortest paths of differing homology. Filled only when
    /// `BallParams::certify_masks`. Expected non-zero on a surface code; recorded, not blocking.
    uint64_t ambiguous_mask_pairs{0};

    uint64_t bytes_targets{0};
    uint64_t bytes_weights{0};
    uint64_t bytes_masks{0};
    uint64_t bytes_paths{0};
    uint64_t bytes_words{0};
    uint64_t bytes_boundary{0};
    uint64_t bytes_offsets{0};
    uint64_t bytes_total{0};

    double compile_wall_seconds{0};
};

/// The offline ball tables of §M2.2: for every node `v` of the detector graph `G`, every node `u`
/// within `R` of it, the **exact** integer distance `d_G(v, u)` in `G`'s discretised metric, and the
/// observable set of the canonical shortest path.
///
/// Storage is CSR pools with no translation-invariance compression: the design deliberately assumes
/// the full table is resident, because the hardware argument is local memory per node.
///
/// **Canonical path.** Neighbours are relaxed in ascending node id, and on an equal tentative
/// distance the predecessor with the smaller node id wins. Every consumer must agree with this, so
/// there is exactly one path per (v, u) and the observable masks are reproducible.
struct BallTables {
    BallParams params;

    /// `MatchingGraph::normalising_constant` of the graph the tables were compiled from. The unit
    /// rule of §0: this is the one constant that converts weights to time units, and the decoder
    /// debug-asserts that it holds the same one.
    double normalising_constant{0};
    /// `to_time_units(R)` and `to_time_units(T_max)`.
    horizon_int r_int{0};
    horizon_int t_max_int{0};
    /// `to_time_units(shell_width)`; 0 when there is a single shell.
    horizon_int shell_width_int{0};
    size_t num_shells{1};

    size_t num_nodes{0};
    size_t num_observables{0};
    /// Hash of the graph the tables were compiled against. Load-time mismatch is a hard error.
    uint64_t graph_hash{0};

    /// Entries of `v`'s ball, `[ball_offsets[v], ball_offsets[v + 1])`, sorted by `(w_int, target)`.
    /// `v` itself is not an entry.
    std::vector<uint64_t> ball_offsets;
    /// Shell boundaries within each ball, stride `num_shells + 1`, absolute indices into
    /// `ball_target`. `ball_shell_offsets[v * (num_shells + 1) + s]` starts shell `s` of `v`.
    std::vector<uint64_t> ball_shell_offsets;
    std::vector<uint32_t> ball_target;
    std::vector<pm::weight_int> ball_w_int;

    /// Observable ids crossed by the canonical path of entry `e`, sorted ascending:
    /// `[ball_mask_offsets[e], ball_mask_offsets[e + 1])` of `ball_mask_ids`.
    std::vector<uint64_t> ball_mask_offsets;
    std::vector<uint32_t> ball_mask_ids;

    /// The canonical path of entry `e` as the node sequence `v = n_0, n_1, ..., n_k = target`; an
    /// *edge* of the path is a consecutive pair. Populated only when `params.store_paths`.
    std::vector<uint64_t> ball_path_offsets;
    std::vector<uint32_t> ball_path_nodes;

    /// `bcost(v)`: the cheapest path from `v` to the boundary. Exact whenever it is `<= R`;
    /// `has_bcost[v]` is 0 when there is no boundary path that short, in which case `H` never gives
    /// `v` a boundary edge (correct, since `H` only uses `bcost` when it is `<= T <= R / 2`).
    std::vector<uint8_t> has_bcost;
    std::vector<pm::weight_int> bcost_w_int;
    std::vector<uint64_t> bcost_mask_offsets;
    std::vector<uint32_t> bcost_mask_ids;
    /// Node sequence from `v` to the node carrying the boundary half-edge. Only when
    /// `params.store_paths`.
    std::vector<uint64_t> bcost_path_offsets;
    std::vector<uint32_t> bcost_path_nodes;

    /// Bitset view: `v`'s ball as a bit per detector, covering the detector index window
    /// `[64 * ball_word_base[v], 64 * (ball_word_base[v] + word_len(v)))`. The shot-time
    /// intersection is then an AND against the syndrome bitset. Detector ids within a ball are not
    /// contiguous, so there is slack; `BallStats::mean_ball_word_len` records how much.
    std::vector<uint64_t> ball_word_offsets;
    std::vector<uint32_t> ball_word_base;
    std::vector<uint64_t> ball_words;
    /// Set bits in all *earlier* words of the same node's window. Turns a hit bit into the rank of
    /// its entry, which is how the BITSET path recovers a weight without a second sorted index.
    std::vector<uint32_t> ball_word_rank;
    /// The ball's entry indices in ascending *target* order — the order the bitset yields hits in —
    /// so that `ball_entry_by_rank[ball_offsets[v] + rank]` is the canonical entry index of the
    /// hit. The pools themselves are in `(w_int, target)` order, which is what the SCAN path wants;
    /// this indirection is what lets both build modes name the *same* entry, and hence the same
    /// weight and the same observable set.
    std::vector<uint64_t> ball_entry_by_rank;

    BallStats stats;

    inline uint64_t ball_begin(size_t v) const {
        return ball_offsets[v];
    }
    inline uint64_t ball_end(size_t v) const {
        return ball_offsets[v + 1];
    }
    inline uint64_t ball_size(size_t v) const {
        return ball_offsets[v + 1] - ball_offsets[v];
    }
    /// Entries of shells `0..s` of `v`, i.e. everything within `(s + 1) * shell_width`.
    inline uint64_t ball_end_through_shell(size_t v, size_t s) const {
        size_t stride = num_shells + 1;
        return ball_shell_offsets[v * stride + std::min(s + 1, num_shells)];
    }
    inline uint64_t word_len(size_t v) const {
        return ball_word_offsets[v + 1] - ball_word_offsets[v];
    }

    /// The horizon, in time units, that these tables support. Decoding above it is a hard error.
    inline horizon_int max_horizon() const {
        return t_max_int;
    }
};

/// FNV-1a over the graph's topology, weights and observable masks. Two graphs with the same hash
/// have the same integer metric, which is exactly what the tables depend on.
uint64_t hash_matching_graph(const pm::MatchingGraph& graph);

/// Compiles the ball tables for `mwpm.flooder.graph` (§M2.2). Embarrassingly parallel over source
/// nodes; `num_threads == 0` means "use the hardware concurrency". The result is independent of the
/// thread count — every source node's Dijkstra is independent and the CSR pools are concatenated in
/// ascending node id.
///
/// When the graph has more than 64 observables the `obs_int` masks are unusable, so `mwpm` must
/// carry a search graph (it stores per-edge observable *id* lists); otherwise this throws.
BallTables compile_ball_tables(const pm::Mwpm& mwpm, const BallParams& params, size_t num_threads = 0);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_MANIFOLD_BALL_TABLES_H
