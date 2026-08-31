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
#include <random>

#include "gtest/gtest.h"

#include "pyrematching/spec_matching/truncation/harvest.h"
#include "pyrematching/spec_matching/truncation/truncated_timeline.h"
#include "tests/spec_matching/spec_matching_test_util.h"

using namespace pm::spec_matching;
using namespace pm::spec_matching::test;

namespace {

/// A random connected graph with a boundary, built directly as a `MatchingGraph` so the weights
/// are exact integers. Weights are even, as `iter_discretized_edges` guarantees for real graphs.
pm::MatchingGraph random_graph(std::mt19937_64& rng, size_t num_nodes) {
    pm::MatchingGraph graph(num_nodes, 1);
    std::uniform_int_distribution<int> weight_dist(1, 12);
    // A spanning path keeps the graph connected, so an odd syndrome is always matchable.
    for (size_t i = 1; i < num_nodes; i++)
        graph.add_edge(i - 1, i, 2 * weight_dist(rng), {});
    // Plus some chords, skipping duplicates (parallel edges are not allowed).
    std::uniform_int_distribution<size_t> node_dist(0, num_nodes - 1);
    std::set<std::pair<size_t, size_t>> present;
    for (size_t i = 1; i < num_nodes; i++)
        present.insert({i - 1, i});
    for (size_t i = 0; i < num_nodes; i++) {
        size_t u = node_dist(rng);
        size_t v = node_dist(rng);
        if (u == v)
            continue;
        auto key = std::minmax(u, v);
        if (!present.insert({key.first, key.second}).second)
            continue;
        graph.add_edge(key.first, key.second, 2 * weight_dist(rng), {});
    }
    graph.add_boundary_edge(0, 2 * weight_dist(rng), {});
    graph.add_boundary_edge(num_nodes - 1, 2 * weight_dist(rng), {});
    return graph;
}

}  // namespace

// Test 11. Random horizons on random syndromes over random graphs. Every harvest invariant has to
// hold in every case, and the instance has to stay reusable throughout.
TEST(SpecMatchingFuzz, RandomHorizonsPreserveHarvestInvariants) {
    std::mt19937_64 rng(2026);
    Harvester harvester;
    std::vector<TightPairRecord> tight_pairs;
    harvester.tight_pairs_out = &tight_pairs;

    for (int trial = 0; trial < 200; trial++) {
        size_t num_nodes = 4 + (size_t)(rng() % 24);
        auto graph = random_graph(rng, num_nodes);
        pm::Mwpm mwpm(pm::GraphFlooder(pm::MatchingGraph(std::move(graph))));

        // The biggest possible dual: every region can grow at most to the graph's total weight.
        pm::cumulative_time_int max_weight = 0;
        for (const auto& node : mwpm.flooder.graph.nodes) {
            for (auto weight : node.neighbor_weights)
                max_weight += weight;
        }

        for (int shot_index = 0; shot_index < 8; shot_index++) {
            std::vector<uint64_t> shot;
            for (size_t node = 0; node < num_nodes; node++) {
                if (rng() % 3 == 0)
                    shot.push_back(node);
            }
            pm::cumulative_time_int horizon = (pm::cumulative_time_int)(rng() % (uint64_t)(max_weight + 1));

            tight_pairs.clear();
            TimelineStatus status;
            HarvestResult harvest;
            ASSERT_NO_THROW({
                status = process_timeline_until_horizon(mwpm, shot, horizon);
                harvest = harvester.harvest_to_obs(mwpm, shot);
            }) << "trial "
               << trial << " shot " << shot_index;

            ASSERT_EQ(mwpm.flooder.horizon, pm::NO_HORIZON);
            ASSERT_TRUE(mwpm.flooder.queue.empty());
            ASSERT_EQ(harvest.residual.size(), (size_t)harvest.num_trees);
            ASSERT_EQ(status == TimelineStatus::TRUNCATED, !harvest.residual.empty());
            ASSERT_TRUE(std::is_sorted(harvest.residual.begin(), harvest.residual.end()));
            ASSERT_EQ(std::adjacent_find(harvest.residual.begin(), harvest.residual.end()), harvest.residual.end());

            for (auto dual_sum : harvest.residual_dual_sum)
                ASSERT_EQ(dual_sum, horizon) << "separation invariant broken in trial " << trial;

            pm::total_weight_int residual_duals = 0;
            for (auto dual_sum : harvest.residual_dual_sum)
                residual_duals += dual_sum;
            ASSERT_EQ(harvest.committed.weight + residual_duals, harvest.dual_sum_at_truncation);

            for (auto residual_defect : harvest.residual)
                ASSERT_TRUE(std::binary_search(shot.begin(), shot.end(), residual_defect));

            for (const auto& pair : tight_pairs) {
                auto distances = dijkstra(mwpm.flooder.graph, pair.inner_node);
                ASSERT_EQ(distances[pair.outer_node], pair.inner_dual_sum + pair.outer_dual_sum)
                    << "a tree-committed pair was not tight, in trial " << trial;
            }
        }
    }
}

// The instance must survive an arbitrary interleaving of truncated and exact decodes, and give the
// same answers as instances that only ever did one or the other.
TEST(SpecMatchingFuzz, InterleavedTruncatedAndExactDecodesAgreeWithFreshInstances) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.01, 120, 97);
    auto mwpm = corpus.to_mwpm();
    auto fresh_mwpm = corpus.to_mwpm();
    Harvester harvester;
    std::mt19937_64 rng(4242);
    auto unit = median_edge_weight(mwpm.flooder.graph);

    for (const auto& shot : corpus.shots) {
        if (rng() % 2 == 0) {
            pm::cumulative_time_int horizon = (pm::cumulative_time_int)(rng() % (uint64_t)(unit * 10 + 1));
            ASSERT_NO_THROW(phase1_decode_to_obs(mwpm, shot, horizon, &harvester));
        }
        auto reused = exact_decode(mwpm, shot);
        auto reference = exact_decode(fresh_mwpm, shot);
        ASSERT_EQ(reused.obs_mask, reference.obs_mask);
        ASSERT_EQ(reused.weight, reference.weight);
    }
}
