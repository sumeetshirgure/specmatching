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

#include <set>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/manifold/ball_decoding.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

Corpus small_corpus() {
    return generate_surface_code_corpus(5, 5, 0.005, 60, 7);
}

BallConfig config_for(const pm::MatchingGraph& graph, double t_edges) {
    double unit = edge_weight_units(graph);
    BallConfig config;
    config.T = t_edges * unit;
    config.ball.T_max = 2.0 * unit;
    config.ball.R = 4.0 * unit;
    return config;
}

bool graphs_equal(const BallGraph& a, const BallGraph& b) {
    if (a.h_to_det != b.h_to_det || a.edges.size() != b.edges.size() ||
        a.boundary_edges.size() != b.boundary_edges.size())
        return false;
    for (size_t k = 0; k < a.edges.size(); k++) {
        if (a.edges[k].i != b.edges[k].i || a.edges[k].j != b.edges[k].j || a.edges[k].w_int != b.edges[k].w_int ||
            a.edges[k].entry != b.edges[k].entry)
            return false;
    }
    for (size_t k = 0; k < a.boundary_edges.size(); k++) {
        if (a.boundary_edges[k].i != b.boundary_edges[k].i || a.boundary_edges[k].w_int != b.boundary_edges[k].w_int)
            return false;
    }
    return true;
}

}  // namespace

// B5. The two intersection strategies must produce the *same* `H` on every shot — otherwise the
// choice of strategy would be a choice of answer.
TEST(ManifoldBallGraph, B5ScanMatchesBitset) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);

    BallGraphArena scan_arena;
    BallGraphArena bitset_arena;
    scan_arena.reset_for_graph(decoder.tables.num_nodes);
    bitset_arena.reset_for_graph(decoder.tables.num_nodes);

    std::vector<uint64_t> seeded;
    size_t shots_with_edges = 0;
    for (const auto& shot : corpus.shots) {
        decoder.compute_seeded_detection_events(shot, seeded);
        build_ball_graph(decoder.tables, seeded, decoder.horizon, scan_arena, BallGraphBuildMode::SCAN);
        build_ball_graph(decoder.tables, seeded, decoder.horizon, bitset_arena, BallGraphBuildMode::BITSET);
        ASSERT_TRUE(graphs_equal(scan_arena.graph, bitset_arena.graph)) << "SCAN and BITSET disagree";
        if (!scan_arena.graph.edges.empty())
            shots_with_edges++;
    }
    ASSERT_GT(shots_with_edges, 0u) << "the corpus produced no edges at all; the comparison is vacuous";
}

// Every pair within `2T` is in `H` and nothing beyond it is — the other half of debug invariant 9,
// checked against the table's own metric rather than against the builder's loop.
TEST(ManifoldBallGraph, EdgeSetIsExactlyThePairsWithinTwoT) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 1.5);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
    const BallTables& tables = decoder.tables;
    pm::cumulative_time_int two_t = 2 * decoder.horizon;

    BallGraphArena arena;
    arena.reset_for_graph(tables.num_nodes);
    std::vector<uint64_t> seeded;
    for (const auto& shot : corpus.shots) {
        decoder.compute_seeded_detection_events(shot, seeded);
        build_ball_graph(tables, seeded, decoder.horizon, arena, BallGraphBuildMode::SCAN);
        const BallGraph& h = arena.graph;

        std::set<std::pair<uint32_t, uint32_t>> present;
        for (const auto& edge : h.edges) {
            ASSERT_LT(edge.i, edge.j);
            ASSERT_LE((pm::cumulative_time_int)edge.w_int, two_t);
            ASSERT_TRUE(present.insert({edge.i, edge.j}).second) << "a pair was emitted twice";
        }
        for (uint32_t i = 0; i < h.num_nodes(); i++) {
            for (uint32_t j = i + 1; j < h.num_nodes(); j++) {
                uint64_t a = h.h_to_det[i];
                uint64_t b = h.h_to_det[j];
                pm::weight_int distance = 0;
                bool within = false;
                for (uint64_t e = tables.ball_begin(a); e < tables.ball_end(a); e++) {
                    if (tables.ball_target[e] == b) {
                        distance = tables.ball_w_int[e];
                        within = (pm::cumulative_time_int)distance <= two_t;
                        break;
                    }
                }
                ASSERT_EQ(within, present.count({i, j}) == 1) << "pair (" << a << ", " << b << ") membership is wrong";
            }
        }
        std::set<uint32_t> boundary;
        for (const auto& edge : h.boundary_edges)
            boundary.insert(edge.i);
        for (uint32_t i = 0; i < h.num_nodes(); i++) {
            uint64_t det = h.h_to_det[i];
            bool expected =
                tables.has_bcost[det] && (pm::cumulative_time_int)tables.bcost_w_int[det] <= decoder.horizon;
            ASSERT_EQ(expected, boundary.count(i) == 1) << "boundary edge for detector " << det << " is wrong";
        }
    }
}

// B7. The arena and the `H`-side `Mwpm` stop growing once the largest shot has been seen, and a
// warmed decoder returns exactly what a fresh one does.
TEST(ManifoldBallGraph, B7ArenaReuse) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    BallDecoder warm = BallDecoder::from_mwpm(std::move(mwpm), config);

    std::vector<CommittedPair> pairs;
    std::vector<HarvestResult> first_pass;
    for (const auto& shot : corpus.shots)
        first_pass.push_back(warm.decode_phase1_to_match_edges(shot, pairs));

    uint64_t arena_grows = warm.arena.grow_events;
    uint64_t mwpm_grows = warm.h_mwpm.grow_events;

    auto fresh_mwpm = corpus.to_mwpm();
    BallDecoder fresh = BallDecoder::from_mwpm(std::move(fresh_mwpm), config);
    for (size_t s = 0; s < corpus.shots.size(); s++) {
        HarvestResult again = warm.decode_phase1_to_match_edges(corpus.shots[s], pairs);
        std::vector<CommittedPair> fresh_pairs;
        HarvestResult from_fresh = fresh.decode_phase1_to_match_edges(corpus.shots[s], fresh_pairs);
        ASSERT_EQ(again.residual, first_pass[s].residual) << "shot " << s << " is not reproducible";
        ASSERT_EQ(again.residual, from_fresh.residual) << "a warmed decoder differs from a fresh one on shot " << s;
        ASSERT_EQ(again.dual_sum_at_truncation, from_fresh.dual_sum_at_truncation);
        ASSERT_EQ(again.committed.weight, from_fresh.committed.weight);
        ASSERT_EQ(pairs.size(), fresh_pairs.size());
    }

    // Debug invariant 11: nothing grew on the second pass.
    ASSERT_EQ(warm.arena.grow_events, arena_grows) << "the ball graph arena kept allocating after warmup";
    ASSERT_EQ(warm.h_mwpm.grow_events, mwpm_grows) << "the H-side Mwpm kept allocating after warmup";
}

// B8. Degenerate shots.
TEST(ManifoldBallGraph, B8EmptyShot) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);

    HarvestResult result = decoder.decode_phase1({});
    ASSERT_TRUE(result.residual.empty());
    ASSERT_EQ(result.num_trees, 0);
    ASSERT_EQ(result.committed.weight, 0);
    ASSERT_EQ(result.committed.obs_mask, 0u);
}

TEST(ManifoldBallGraph, B8ZeroHorizonCommitsNothing) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    config.T = 0;
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
    ASSERT_EQ(decoder.horizon, 0);

    std::vector<uint64_t> seeded;
    for (const auto& shot : corpus.shots) {
        HarvestResult result = decoder.decode_phase1(shot);
        decoder.compute_seeded_detection_events(shot, seeded);
        // At `T = 0` nothing can have interacted, so the residual is exactly the seeded events.
        ASSERT_EQ(result.residual, seeded);
        ASSERT_EQ(result.committed.weight, 0);
        ASSERT_EQ(result.committed_pairs_frozen, 0);
        ASSERT_EQ(result.committed_pairs_tree, 0);
        ASSERT_EQ(result.committed_boundary, 0);
    }
}

TEST(ManifoldBallGraph, B8MutuallyDistantDefectsGiveAnEdgelessGraph) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 1.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
    const BallTables& tables = decoder.tables;
    pm::cumulative_time_int two_t = 2 * decoder.horizon;

    // Greedily pick defects that are pairwise beyond `2T` and beyond `T` of the boundary, so that
    // `H` has no edges at all and every one of them has to survive as its own tree.
    std::vector<uint64_t> defects;
    for (uint64_t v = 0; v < tables.num_nodes; v++) {
        if (tables.has_bcost[v] && (pm::cumulative_time_int)tables.bcost_w_int[v] <= decoder.horizon)
            continue;
        bool clear = true;
        for (uint64_t u : defects) {
            for (uint64_t e = tables.ball_begin(v); e < tables.ball_end(v); e++) {
                if (tables.ball_target[e] == u && (pm::cumulative_time_int)tables.ball_w_int[e] <= two_t) {
                    clear = false;
                    break;
                }
            }
            if (!clear)
                break;
        }
        if (clear)
            defects.push_back(v);
    }
    ASSERT_GE(defects.size(), 2u) << "could not build a mutually distant defect set";

    HarvestResult result = decoder.decode_phase1(defects);
    ASSERT_EQ(result.residual, defects);
    ASSERT_EQ(result.num_trees, (int)defects.size());
    ASSERT_EQ(result.committed.weight, 0);
}

// Asking the ball front end for an infinite horizon is a hard error, not a silent truncation: `H`
// would then be the complete defect graph, which no ball radius covers.
TEST(ManifoldBallGraph, InfiniteHorizonIsRejected) {
    Corpus corpus = small_corpus();
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
    BallGraphArena arena;
    arena.reset_for_graph(decoder.tables.num_nodes);
    ASSERT_THROW(
        build_ball_graph(decoder.tables, {0}, pm::NO_HORIZON, arena, BallGraphBuildMode::SCAN), std::invalid_argument);
}
