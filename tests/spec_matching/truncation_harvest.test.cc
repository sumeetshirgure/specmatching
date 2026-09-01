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

#include <algorithm>

#include "gtest/gtest.h"

#include "specmatching/spec_matching/truncation/harvest.h"
#include "specmatching/spec_matching/truncation/truncated_timeline.h"
#include "tests/spec_matching/spec_matching_test_util.h"

using namespace pm::spec_matching;
using namespace pm::spec_matching::test;

// Test 2, part 1. The separation invariant: a residual defect's dual is exactly the horizon.
//
// This is the property everything downstream rests on. If it ever fires the Phase-2 error bound is
// void, and that is a release blocker rather than a flaky test.
TEST(SpecMatchingHarvest, ResidualDefectsHaveDualSumExactlyEqualToTheHorizon) {
    std::vector<Corpus> corpora;
    corpora.push_back(load_surface_code_d13(60));
    corpora.push_back(load_surface_code_d13_negative_weights(60));
    corpora.push_back(load_toric_code_d5(60));

    for (const Corpus& corpus : corpora) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        size_t residual_defects_seen = 0;
        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3, 5, 8}) {
                pm::cumulative_time_int horizon = unit * multiple;
                auto decoded = phase1_decode_to_obs(mwpm, shot, horizon, &harvester);
                for (size_t i = 0; i < decoded.harvest.residual.size(); i++) {
                    ASSERT_EQ(decoded.harvest.residual_dual_sum[i], horizon)
                        << corpus.name << ": Y(u) != T for residual defect " << decoded.harvest.residual[i]
                        << " at T=" << horizon;
                    residual_defects_seen++;
                }
            }
        }
        ASSERT_GT(residual_defects_seen, 0u) << corpus.name << ": no truncation happened, nothing was tested";
    }
}

// Test 2, part 2. On a graph small enough for an exact Dijkstra, residual defects are pairwise at
// least 2T apart and at least T from the boundary — the geometric form of the same invariant.
TEST(SpecMatchingHarvest, ResidualDefectsAreMutuallySeparatedByTwiceTheHorizon) {
    Corpus corpus = generate_surface_code_corpus(9, 9, 0.01, 100, 7);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    size_t pairs_checked = 0;
    for (const auto& shot : corpus.shots) {
        for (int multiple : {1, 2, 4}) {
            pm::cumulative_time_int horizon = unit * multiple;
            auto decoded = phase1_decode_to_obs(mwpm, shot, horizon, &harvester);
            const auto& residual = decoded.harvest.residual;
            if (residual.empty())
                continue;

            for (size_t i = 0; i < residual.size(); i++) {
                auto distances = dijkstra(mwpm.flooder.graph, residual[i]);
                for (size_t j = i + 1; j < residual.size(); j++) {
                    ASSERT_GE(distances[residual[j]], 2 * horizon)
                        << "residual defects " << residual[i] << " and " << residual[j] << " are closer than 2T";
                    pairs_checked++;
                }
                auto to_boundary = boundary_cost(mwpm.flooder.graph, residual[i]);
                if (to_boundary >= 0)
                    ASSERT_GE(to_boundary, horizon)
                        << "residual defect " << residual[i] << " is closer than T to the boundary";
            }
        }
    }
    ASSERT_GT(pairs_checked, 0u) << "no residual pair was ever produced, so nothing was tested";
}

// Invariant 3. Harvest classifies every detection event exactly once: committed and residual
// partition the seeded detection events, and there is exactly one residual defect per tree.
TEST(SpecMatchingHarvest, CommittedAndResidualPartitionTheDetectionEvents) {
    Corpus corpus = load_surface_code_d13(40);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    for (const auto& shot : corpus.shots) {
        for (int multiple : {0, 1, 3, 6}) {
            pm::cumulative_time_int horizon = unit * multiple;
            auto seeded = expected_seeded_detection_events(mwpm, shot);

            auto status = process_timeline_until_horizon(mwpm, shot, horizon);
            std::vector<pm::CompressedEdge> match_edges;
            auto harvest = harvester.harvest_to_match_edges(mwpm, shot, match_edges);

            ASSERT_EQ(harvest.residual.size(), (size_t)harvest.num_trees);
            ASSERT_EQ(status == TimelineStatus::TRUNCATED, !harvest.residual.empty());

            std::vector<uint64_t> classified = harvest.residual;
            for (const auto& edge : match_edges) {
                classified.push_back((uint64_t)(edge.loc_from - mwpm.flooder.graph.nodes.data()));
                if (edge.loc_to != nullptr)
                    classified.push_back((uint64_t)(edge.loc_to - mwpm.flooder.graph.nodes.data()));
            }
            std::sort(classified.begin(), classified.end());
            ASSERT_EQ(std::adjacent_find(classified.begin(), classified.end()), classified.end())
                << "a detection event was classified twice";
            ASSERT_EQ(classified, seeded) << "committed + residual is not a partition of the seeded events";
        }
    }
}

// Test 7. The truncated solution can only be worse than the exact one: committing early and then
// solving the residual exactly cannot beat solving the whole thing exactly.
TEST(SpecMatchingHarvest, CommittedPlusExactResidualUpperBoundsTheExactOptimum) {
    Corpus corpus = generate_surface_code_corpus(7, 7, 0.005, 150, 11);
    auto mwpm = corpus.to_mwpm();
    auto residual_mwpm = corpus.to_mwpm();
    auto reference_mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    for (const auto& shot : corpus.shots) {
        auto exact_full = exact_decode(reference_mwpm, shot);
        for (int multiple : {1, 2, 4, 8}) {
            auto decoded = phase1_decode_to_obs(mwpm, shot, unit * multiple, &harvester);
            auto exact_residual = exact_decode(residual_mwpm, decoded.harvest.residual);
            ASSERT_GE(decoded.harvest.committed.weight + exact_residual.weight, exact_full.weight)
                << "the spec-matching split beat the exact optimum, which is impossible";
            // The dual solution at truncation is feasible, so it lower bounds the exact optimum.
            ASSERT_LE(decoded.harvest.dual_sum_at_truncation, exact_full.weight);
        }
    }
}

// Test 9 / invariant 5. Every pair committed out of a surviving tree is joined by a tight edge:
// the graph distance between its endpoints is exactly the sum of their duals.
TEST(SpecMatchingHarvest, TreeCommittedPairsAreTight) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.008, 120, 23);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    std::vector<TightPairRecord> tight_pairs;
    harvester.tight_pairs_out = &tight_pairs;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    size_t pairs_checked = 0;
    for (const auto& shot : corpus.shots) {
        for (int multiple : {1, 2, 4}) {
            tight_pairs.clear();
            phase1_decode_to_obs(mwpm, shot, unit * multiple, &harvester);
            for (const auto& pair : tight_pairs) {
                auto distances = dijkstra(mwpm.flooder.graph, pair.inner_node);
                ASSERT_EQ(distances[pair.outer_node], pair.inner_dual_sum + pair.outer_dual_sum)
                    << "tree-committed pair (" << pair.inner_node << ", " << pair.outer_node << ") is not tight";
                pairs_checked++;
            }
        }
    }
    ASSERT_GT(pairs_checked, 0u) << "no tree pair was ever committed, so nothing was tested";
}

// The dual solution at truncation is feasible, so it lower bounds the exact optimum: this is the
// per-shot certificate the M3 accuracy sweep reports. It also splits exactly, which is the
// accounting statement of the commit policy: every region's dual either paid for a committed pair
// or belongs to a residual defect, and nothing is counted twice or dropped.
TEST(SpecMatchingHarvest, DualSumSplitsExactlyBetweenCommittedWeightAndResidualDefects) {
    Corpus corpus = load_surface_code_d13(40);
    auto mwpm = corpus.to_mwpm();
    auto reference_mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    for (const auto& shot : corpus.shots) {
        auto exact_full = exact_decode(reference_mwpm, shot);
        for (int multiple : {1, 3, 6, 12}) {
            auto decoded = phase1_decode_to_obs(mwpm, shot, unit * multiple, &harvester);
            pm::total_weight_int residual_duals = 0;
            for (auto dual_sum : decoded.harvest.residual_dual_sum)
                residual_duals += dual_sum;
            ASSERT_EQ(decoded.harvest.committed.weight + residual_duals, decoded.harvest.dual_sum_at_truncation);
            ASSERT_LE(decoded.harvest.dual_sum_at_truncation, exact_full.weight);
        }
    }
}
