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

#include <algorithm>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

std::vector<Corpus> identity_corpora() {
    std::vector<Corpus> corpora;
    corpora.push_back(load_surface_code_d13(60));
    corpora.push_back(load_surface_code_d13_negative_weights(60));
    corpora.push_back(load_toric_code_d5(60));
    return corpora;
}

std::vector<std::pair<int64_t, int64_t>> to_sorted_index_pairs(
    const std::vector<pm::CompressedEdge>& match_edges, const pm::MatchingGraph& graph) {
    std::vector<std::pair<int64_t, int64_t>> pairs;
    pairs.reserve(match_edges.size());
    for (const auto& edge : match_edges) {
        int64_t from = edge.loc_from - graph.nodes.data();
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)(edge.loc_to - graph.nodes.data());
        pairs.push_back(from <= to ? std::make_pair(from, to) : std::make_pair(to, from));
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

}  // namespace

// Test 1. The truncated timeline at T = infinity is the stock decoder, bit for bit. This is the
// gate that says the horizon machinery changed nothing on the path everybody else uses.
TEST(TwoPhaseTruncation, InfiniteHorizonIsBitExactObsFlavor) {
    for (const Corpus& corpus : identity_corpora()) {
        auto mwpm = corpus.to_mwpm();
        auto reference_mwpm = corpus.to_mwpm();
        Harvester harvester;

        for (const auto& shot : corpus.shots) {
            auto expected = exact_decode(reference_mwpm, shot);
            auto actual = phase1_decode_to_obs(mwpm, shot, pm::NO_HORIZON, &harvester);

            ASSERT_EQ(actual.status, TimelineStatus::COMPLETE) << corpus.name;
            ASSERT_TRUE(actual.harvest.residual.empty()) << corpus.name;
            ASSERT_EQ(actual.harvest.num_trees, 0) << corpus.name;
            ASSERT_EQ(actual.obs_mask, expected.obs_mask) << corpus.name;
            ASSERT_EQ(actual.weight, expected.weight) << corpus.name;
            // A completed timeline's dual solution is tight against its own primal.
            ASSERT_EQ(actual.harvest.dual_sum_at_truncation, actual.harvest.committed.weight) << corpus.name;
        }
    }
}

TEST(TwoPhaseTruncation, InfiniteHorizonIsBitExactMatchEdgesFlavor) {
    // The negative-weight corpus is excluded: `decode_detection_events_to_match_edges` refuses
    // graphs with negative weight edges, so there is nothing to compare against.
    std::vector<Corpus> corpora;
    corpora.push_back(load_surface_code_d13(40));
    corpora.push_back(load_toric_code_d5(40));
    for (const Corpus& corpus : corpora) {
        auto mwpm = corpus.to_mwpm();
        auto reference_mwpm = corpus.to_mwpm();
        Harvester harvester;

        for (const auto& shot : corpus.shots) {
            pm::decode_detection_events_to_match_edges(reference_mwpm, shot);
            auto expected = to_sorted_index_pairs(reference_mwpm.flooder.match_edges, reference_mwpm.flooder.graph);
            reference_mwpm.reset();

            auto status = process_timeline_until_horizon(mwpm, shot, pm::NO_HORIZON);
            std::vector<pm::CompressedEdge> match_edges;
            auto harvest = harvester.harvest_to_match_edges(mwpm, shot, match_edges);
            auto actual = to_sorted_index_pairs(match_edges, mwpm.flooder.graph);

            ASSERT_EQ(status, TimelineStatus::COMPLETE) << corpus.name;
            ASSERT_TRUE(harvest.residual.empty()) << corpus.name;
            ASSERT_EQ(actual, expected) << corpus.name;
        }
    }
}

// Test 4. At T = 0 nothing has had time to happen: every seeded detection event stays exposed.
TEST(TwoPhaseTruncation, ZeroHorizonCommitsNothing) {
    Corpus corpus = load_surface_code_d13(30);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;

    for (const auto& shot : corpus.shots) {
        auto expected_residual = expected_seeded_detection_events(mwpm, shot);
        auto decoded = phase1_decode_to_obs(mwpm, shot, 0, &harvester);

        if (expected_residual.empty()) {
            ASSERT_EQ(decoded.status, TimelineStatus::COMPLETE);
        } else {
            ASSERT_EQ(decoded.status, TimelineStatus::TRUNCATED);
        }
        ASSERT_EQ(decoded.harvest.residual, expected_residual);
        ASSERT_EQ(decoded.harvest.committed.weight, 0);
        ASSERT_EQ(decoded.harvest.committed_pairs_frozen, 0);
        ASSERT_EQ(decoded.harvest.committed_pairs_tree, 0);
        ASSERT_EQ(decoded.harvest.committed_boundary, 0);
        ASSERT_EQ(decoded.harvest.dual_sum_at_truncation, 0);
        for (auto dual_sum : decoded.harvest.residual_dual_sum)
            ASSERT_EQ(dual_sum, 0);
    }
}

// Test 5. Growing the horizon can only commit more, so the residual can only shrink.
TEST(TwoPhaseTruncation, ResidualSizeIsMonotonicInHorizon) {
    Corpus corpus = load_surface_code_d13(30);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    std::vector<pm::cumulative_time_int> horizons;
    for (int multiple = 0; multiple <= 12; multiple++)
        horizons.push_back(unit * multiple);

    for (const auto& shot : corpus.shots) {
        size_t previous_residual_size = SIZE_MAX;
        for (auto horizon : horizons) {
            auto decoded = phase1_decode_to_obs(mwpm, shot, horizon, &harvester);
            ASSERT_LE(decoded.harvest.residual.size(), previous_residual_size)
                << "residual grew when the horizon grew, at T=" << horizon;
            previous_residual_size = decoded.harvest.residual.size();
        }
    }
}

// Test 6. Truncation is a normal outcome, not an error, and it leaves the instance reusable.
TEST(TwoPhaseTruncation, TruncationNeverThrowsAndLeavesInstanceReusable) {
    Corpus corpus = load_surface_code_d13(40);
    auto mwpm = corpus.to_mwpm();
    auto fresh_mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    size_t truncated_shots = 0;
    for (size_t i = 0; i + 1 < corpus.shots.size(); i++) {
        // Truncate one shot...
        Phase1Decode truncated;
        ASSERT_NO_THROW({ truncated = phase1_decode_to_obs(mwpm, corpus.shots[i], unit * 2, &harvester); });
        if (truncated.status == TimelineStatus::TRUNCATED)
            truncated_shots++;
        ASSERT_EQ(truncated.status == TimelineStatus::TRUNCATED, !truncated.harvest.residual.empty());
        ASSERT_EQ(mwpm.flooder.horizon, pm::NO_HORIZON) << "the sentinel horizon was not restored";

        // ...then decode the next one exactly, and demand the same answer a pristine instance gives.
        auto reused = exact_decode(mwpm, corpus.shots[i + 1]);
        auto reference = exact_decode(fresh_mwpm, corpus.shots[i + 1]);
        ASSERT_EQ(reused.obs_mask, reference.obs_mask) << "shot " << i + 1 << " after a truncated shot";
        ASSERT_EQ(reused.weight, reference.weight) << "shot " << i + 1 << " after a truncated shot";
    }
    ASSERT_GT(truncated_shots, 0u) << "the horizon was too generous for this test to test anything";
}

// Truncating and then harvesting must leave the flooder queue empty, or the next shot's stock
// preamble refuses to start.
TEST(TwoPhaseTruncation, HarvestLeavesTheQueueEmpty) {
    Corpus corpus = load_surface_code_d13(20);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    for (const auto& shot : corpus.shots) {
        phase1_decode_to_obs(mwpm, shot, unit * 3, &harvester);
        ASSERT_TRUE(mwpm.flooder.queue.empty());
        ASSERT_EQ(mwpm.node_arena.allocated.size(), mwpm.node_arena.available.size());
    }
}

// M1 exit checkpoint: a debug counter confirms that nothing past the horizon is ever enqueued.
// The counters only exist in debug builds, which is where the gtest target is built.
TEST(TwoPhaseTruncation, NothingPastTheHorizonIsEverEnqueued) {
#ifdef NDEBUG
    GTEST_SKIP() << "horizon gate counters are compiled out of release builds";
#else
    Corpus corpus = load_surface_code_d13(30);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    uint64_t total_rejected = 0;
    for (const auto& shot : corpus.shots) {
        for (auto horizon : {unit, unit * 3, unit * 8}) {
            pm::horizon_gate_stats.clear();
            phase1_decode_to_obs(mwpm, shot, horizon, &harvester);
            ASSERT_LE(pm::horizon_gate_stats.max_admitted_time, horizon)
                << "an event past the horizon reached the queue";
            total_rejected += pm::horizon_gate_stats.rejected;
        }
    }
    ASSERT_GT(total_rejected, 0u) << "the gate never fired, so this test proved nothing";
    pm::horizon_gate_stats.clear();
#endif
}

// Invariant 1, from the other side: with the sentinel horizon the gate must never reject anything.
TEST(TwoPhaseTruncation, SentinelHorizonRejectsNothing) {
#ifdef NDEBUG
    GTEST_SKIP() << "horizon gate counters are compiled out of release builds";
#else
    Corpus corpus = load_surface_code_d13(20);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;

    pm::horizon_gate_stats.clear();
    for (const auto& shot : corpus.shots)
        phase1_decode_to_obs(mwpm, shot, pm::NO_HORIZON, &harvester);
    ASSERT_EQ(pm::horizon_gate_stats.rejected, 0u);
    pm::horizon_gate_stats.clear();
#endif
}
