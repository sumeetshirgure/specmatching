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

#include <cstdio>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/manifold/ball_serialize.h"
#include "pyrematching/two_phase/manifold/ball_tables.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

/// A graph small enough that an all-pairs sweep and a per-entry path check are affordable, which is
/// what B1 and B3 need.
Corpus small_corpus() {
    return generate_surface_code_corpus(5, 5, 0.005, 40, 2026);
}

BallParams params_for(const pm::MatchingGraph& graph, double t_max_edges, double r_edges, bool paths, bool certify) {
    double unit = edge_weight_units(graph);
    BallParams params;
    params.T_max = t_max_edges * unit;
    params.R = r_edges * unit;
    params.store_paths = paths;
    params.certify_masks = certify;
    return params;
}

bool tables_equal(const BallTables& a, const BallTables& b) {
    return a.params == b.params && a.normalising_constant == b.normalising_constant && a.r_int == b.r_int &&
           a.t_max_int == b.t_max_int && a.shell_width_int == b.shell_width_int && a.num_shells == b.num_shells &&
           a.num_nodes == b.num_nodes && a.num_observables == b.num_observables && a.graph_hash == b.graph_hash &&
           a.ball_offsets == b.ball_offsets && a.ball_shell_offsets == b.ball_shell_offsets &&
           a.ball_target == b.ball_target && a.ball_w_int == b.ball_w_int &&
           a.ball_mask_offsets == b.ball_mask_offsets && a.ball_mask_ids == b.ball_mask_ids &&
           a.ball_path_offsets == b.ball_path_offsets && a.ball_path_nodes == b.ball_path_nodes &&
           a.has_bcost == b.has_bcost && a.bcost_w_int == b.bcost_w_int &&
           a.bcost_mask_offsets == b.bcost_mask_offsets && a.bcost_mask_ids == b.bcost_mask_ids &&
           a.bcost_path_offsets == b.bcost_path_offsets && a.bcost_path_nodes == b.bcost_path_nodes &&
           a.ball_word_offsets == b.ball_word_offsets && a.ball_word_base == b.ball_word_base &&
           a.ball_words == b.ball_words && a.ball_word_rank == b.ball_word_rank &&
           a.ball_entry_by_rank == b.ball_entry_by_rank;
}

}  // namespace

// B1. Every stored `w_int` is the integer-metric shortest path distance in `G`, and every node
// within `R` is in the ball. Checked against an independent, unbounded Dijkstra — the ball
// compiler's own bounded sweep is not allowed to be its own oracle.
TEST(ManifoldBallTables, B1TableCorrectness) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    const auto& graph = mwpm.flooder.graph;
    BallTables tables = compile_ball_tables(mwpm, params_for(graph, 2.0, 4.0, false, false));

    ASSERT_GT(tables.stats.total_entries, 0u);
    for (size_t v = 0; v < graph.nodes.size(); v++) {
        auto reference = dijkstra(graph, v);

        // Everything stored is right...
        for (uint64_t e = tables.ball_begin(v); e < tables.ball_end(v); e++) {
            uint32_t u = tables.ball_target[e];
            ASSERT_NE(u, v);
            ASSERT_EQ((pm::total_weight_int)tables.ball_w_int[e], reference[u]) << "v=" << v << " u=" << u;
        }
        // ...and nothing within `R` is missing.
        size_t within = 0;
        for (size_t u = 0; u < graph.nodes.size(); u++) {
            if (u != v && reference[u] <= tables.r_int)
                within++;
        }
        ASSERT_EQ(within, tables.ball_size(v)) << "ball at node " << v << " has the wrong membership";

        // `bcost` is exact whenever it is within `R`, which is all `H` ever asks of it.
        pm::total_weight_int expected_bcost = boundary_cost(graph, v);
        if (expected_bcost >= 0 && expected_bcost <= tables.r_int) {
            ASSERT_TRUE(tables.has_bcost[v]) << "node " << v << " lost its boundary edge";
            ASSERT_EQ((pm::total_weight_int)tables.bcost_w_int[v], expected_bcost) << "node " << v;
        } else {
            ASSERT_FALSE(tables.has_bcost[v]) << "node " << v << " gained a boundary edge beyond R";
        }
    }
}

// B2. Symmetry, ordering and shells. A violation of symmetry means `H` would depend on which
// endpoint of a pair the shot-time intersection walked from.
TEST(ManifoldBallTables, B2SymmetryAndShells) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    const auto& graph = mwpm.flooder.graph;
    double unit = edge_weight_units(graph);
    BallParams params = params_for(graph, 2.0, 4.0, false, false);
    params.shell_width = unit;  // four shells across the ball
    BallTables tables = compile_ball_tables(mwpm, params);
    ASSERT_GT(tables.num_shells, 1u);

    size_t stride = tables.num_shells + 1;
    for (size_t v = 0; v < graph.nodes.size(); v++) {
        // Sorted by `(w_int, target)`, which is what lets the SCAN path stop at the first entry
        // past `2T`.
        for (uint64_t e = tables.ball_begin(v) + 1; e < tables.ball_end(v); e++) {
            bool ordered =
                tables.ball_w_int[e - 1] < tables.ball_w_int[e] ||
                (tables.ball_w_int[e - 1] == tables.ball_w_int[e] && tables.ball_target[e - 1] < tables.ball_target[e]);
            ASSERT_TRUE(ordered) << "ball entries out of order at node " << v;
        }

        // Membership and weight symmetry.
        for (uint64_t e = tables.ball_begin(v); e < tables.ball_end(v); e++) {
            uint32_t u = tables.ball_target[e];
            bool found = false;
            for (uint64_t f = tables.ball_begin(u); f < tables.ball_end(u); f++) {
                if (tables.ball_target[f] == v) {
                    ASSERT_EQ(tables.ball_w_int[f], tables.ball_w_int[e]);
                    found = true;
                }
            }
            ASSERT_TRUE(found) << u << " is in B(" << v << ") but not the other way round";
        }

        // Shells partition the ball.
        ASSERT_EQ(tables.ball_shell_offsets[v * stride], tables.ball_begin(v));
        ASSERT_EQ(tables.ball_shell_offsets[v * stride + tables.num_shells], tables.ball_end(v));
        for (size_t s = 0; s < tables.num_shells; s++) {
            uint64_t begin = tables.ball_shell_offsets[v * stride + s];
            uint64_t end = tables.ball_shell_offsets[v * stride + s + 1];
            ASSERT_LE(begin, end);
            for (uint64_t e = begin; e < end; e++) {
                pm::cumulative_time_int w = tables.ball_w_int[e];
                ASSERT_GE(w, (pm::cumulative_time_int)s * tables.shell_width_int);
                ASSERT_LT(w, (pm::cumulative_time_int)(s + 1) * tables.shell_width_int);
            }
        }
    }
}

// B3. Every stored mask is the XOR of its stored path's edge observables, recomputed
// independently, and every stored path is a valid `G` path of exactly the stored weight.
TEST(ManifoldBallTables, B3MaskConsistency) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    const auto& graph = mwpm.flooder.graph;
    BallTables tables = compile_ball_tables(mwpm, params_for(graph, 2.0, 4.0, true, true));
    ASSERT_LE(graph.num_observables, 64u);

    auto mask_of_slice = [&](const std::vector<uint32_t>& ids, uint64_t begin, uint64_t end) {
        pm::obs_int mask = 0;
        for (uint64_t k = begin; k < end; k++)
            mask ^= (pm::obs_int)1 << ids[k];
        return mask;
    };

    for (size_t v = 0; v < graph.nodes.size(); v++) {
        for (uint64_t e = tables.ball_begin(v); e < tables.ball_end(v); e++) {
            uint64_t path_begin = tables.ball_path_offsets[e];
            uint64_t path_end = tables.ball_path_offsets[e + 1];
            ASSERT_GE(path_end - path_begin, 2u) << "a path must have at least a source and a target";
            ASSERT_EQ(tables.ball_path_nodes[path_begin], v);
            ASSERT_EQ(tables.ball_path_nodes[path_end - 1], tables.ball_target[e]);

            pm::total_weight_int weight = 0;
            pm::obs_int mask = 0;
            for (uint64_t k = path_begin; k + 1 < path_end; k++) {
                uint32_t a = tables.ball_path_nodes[k];
                uint32_t b = tables.ball_path_nodes[k + 1];
                size_t index = graph.nodes[a].index_of_neighbor(const_cast<pm::DetectorNode*>(&graph.nodes[b]));
                ASSERT_LT(index, graph.nodes[a].neighbors.size()) << "stored path uses a non-existent edge";
                weight += graph.nodes[a].neighbor_weights[index];
                mask ^= graph.nodes[a].neighbor_observables[index];
            }
            ASSERT_EQ(weight, (pm::total_weight_int)tables.ball_w_int[e]) << "stored path is not of the stored weight";
            ASSERT_EQ(
                mask,
                mask_of_slice(tables.ball_mask_ids, tables.ball_mask_offsets[e], tables.ball_mask_offsets[e + 1]));
        }

        if (!tables.has_bcost[v])
            continue;
        uint64_t path_begin = tables.bcost_path_offsets[v];
        uint64_t path_end = tables.bcost_path_offsets[v + 1];
        ASSERT_GE(path_end - path_begin, 1u);
        ASSERT_EQ(tables.bcost_path_nodes[path_begin], v);
        pm::total_weight_int weight = 0;
        pm::obs_int mask = 0;
        for (uint64_t k = path_begin; k + 1 < path_end; k++) {
            uint32_t a = tables.bcost_path_nodes[k];
            uint32_t b = tables.bcost_path_nodes[k + 1];
            size_t index = graph.nodes[a].index_of_neighbor(const_cast<pm::DetectorNode*>(&graph.nodes[b]));
            weight += graph.nodes[a].neighbor_weights[index];
            mask ^= graph.nodes[a].neighbor_observables[index];
        }
        // The half-edge into the boundary itself closes the path.
        const auto& last = graph.nodes[tables.bcost_path_nodes[path_end - 1]];
        ASSERT_FALSE(last.neighbors.empty());
        ASSERT_EQ(last.neighbors[0], nullptr) << "a boundary path must end on a node with a boundary half-edge";
        weight += last.neighbor_weights[0];
        mask ^= last.neighbor_observables[0];
        ASSERT_EQ(weight, (pm::total_weight_int)tables.bcost_w_int[v]);
        ASSERT_EQ(
            mask, mask_of_slice(tables.bcost_mask_ids, tables.bcost_mask_offsets[v], tables.bcost_mask_offsets[v + 1]));
    }

    // Reported, never blocking (§M2.2). A non-zero count means some pair inside `R` has two
    // shortest paths of differing homology, which needs a cycle of weight `<= 2R` that flips an
    // observable — impossible once the code distance exceeds `2R`, but perfectly possible on the
    // `d = 5` fixture used here, where `R` is four edge weights. What matters downstream is §M2.6
    // level 2, which checks the observables the decoder actually emits.
    std::printf(
        "[ INFO     ] ambiguous_mask_pairs = %llu of %llu entries (d = 5, R = 4 edge weights)\n",
        (unsigned long long)tables.stats.ambiguous_mask_pairs,
        (unsigned long long)tables.stats.total_entries);
}

// B6. Compilation is deterministic — including across thread counts, which is the only way the
// parallel sweep is allowed to be parallel — and the artifact round-trips byte for byte.
TEST(ManifoldBallTables, B6DeterminismAndRoundTrip) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    const auto& graph = mwpm.flooder.graph;
    BallParams params = params_for(graph, 2.0, 4.0, true, false);
    params.shell_width = edge_weight_units(graph);

    BallTables single = compile_ball_tables(mwpm, params, 1);
    BallTables parallel = compile_ball_tables(mwpm, params, 4);
    ASSERT_TRUE(tables_equal(single, parallel)) << "ball table compilation depends on the thread count";

    std::string path = std::tmpnam(nullptr);
    save_ball_tables(single, path);
    BallTables loaded = load_ball_tables(path, graph, params);
    std::remove(path.c_str());
    ASSERT_TRUE(tables_equal(single, loaded)) << "ball artifact did not round-trip";

    // A mismatched artifact is a hard error, never a best-effort load.
    BallParams other = params;
    other.T_max *= 0.5;
    other.R *= 0.5;
    save_ball_tables(single, path);
    ASSERT_THROW(load_ball_tables(path, graph, other), std::runtime_error);
    std::remove(path.c_str());
}

// `R >= 2 * T_max` is checked in release builds too: an undersized `R` changes the answer with no
// crash to point at it.
TEST(ManifoldBallTables, UndersizedRadiusIsRejected) {
    BallParams params;
    params.T_max = 2.0;
    params.R = 3.9;
    ASSERT_THROW(params.validate(), std::invalid_argument);
    params.R = 4.0;
    ASSERT_NO_THROW(params.validate());
}
