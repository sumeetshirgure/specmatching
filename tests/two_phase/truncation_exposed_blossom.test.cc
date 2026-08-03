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
#include <set>

#include "gtest/gtest.h"

#include "pyrematching/sparse_blossom/matcher/alternating_tree.h"
#include "pyrematching/two_phase/truncation/exposed_blossom.h"
#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

/// The exposed (unmatched) region of every alternating tree that survived truncation.
std::vector<pm::GraphFillRegion*> surviving_exposed_regions(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    std::set<const pm::AltTreeNode*> seen_roots;
    std::vector<pm::GraphFillRegion*> exposed;
    for (uint64_t det : detection_events) {
        if (det >= mwpm.flooder.graph.nodes.size())
            continue;
        auto& node = mwpm.flooder.graph.nodes[det];
        if (node.region_that_arrived == nullptr || node.region_that_arrived_top->alt_tree_node == nullptr)
            continue;
        const pm::AltTreeNode* root = node.region_that_arrived_top->alt_tree_node->find_root();
        if (seen_roots.insert(root).second)
            exposed.push_back(root->outer_region);
    }
    return exposed;
}

/// A three-defect graph engineered so that a single blossom forms at a known time, with a known
/// base, and stays exposed at the horizon.
///
///   0 --4-- 1        0 --8-- 2        1 --6-- 2        0 --40-- boundary
///
/// Timeline: r0 and r1 collide at t=2 and match. r2 reaches the matched r1 at t=4, which turns the
/// pair into an alternating tree below the r2 root (r1 shrinking, r0 growing). r2 and r0 then
/// collide at t=5, which closes the odd cycle {r2, r1, r0} into a blossom at the tree root. The
/// blossom's base is r2: it is the region that has been growing, unmatched, since t=0.
pm::Mwpm build_single_blossom_fixture() {
    pm::MatchingGraph graph(3, 1);
    graph.add_edge(0, 1, 4, {});
    graph.add_edge(1, 2, 6, {});
    graph.add_edge(0, 2, 8, {});
    graph.add_boundary_edge(0, 40, {});
    return pm::Mwpm(pm::GraphFlooder(std::move(graph)));
}

/// The fixture above with a second odd cycle wrapped around it, so that the exposed root region is
/// a blossom whose base is itself a blossom.
///
/// r3 and r4 collide at t=2 and match. The first blossom B = {r2, r1, r0} forms at t=5 exactly as
/// above. B reaches the matched r3 at t=10, making r3 shrink and r4 grow, and then collides with
/// r4 at t=11, closing the cycle {B, r3, r4} into a second, nested blossom. Its base is B, whose
/// base is r2 — so the exposed defect is still detector 2, two levels down.
pm::Mwpm build_nested_blossom_fixture() {
    pm::MatchingGraph graph(5, 1);
    graph.add_edge(0, 1, 4, {});
    graph.add_edge(1, 2, 6, {});
    graph.add_edge(0, 2, 8, {});
    graph.add_edge(2, 3, 12, {});
    graph.add_edge(0, 4, 12, {});
    graph.add_edge(3, 4, 4, {});
    graph.add_boundary_edge(0, 40, {});
    return pm::Mwpm(pm::GraphFlooder(std::move(graph)));
}

}  // namespace

// Test 3, hand-built fixture: a single exposed blossom with a base that is known in advance.
TEST(TwoPhaseExposedBlossom, SingleBlossomExposesItsBase) {
    auto mwpm = build_single_blossom_fixture();
    std::vector<uint64_t> detection_events = {0, 1, 2};
    const pm::cumulative_time_int horizon = 6;

    auto status = process_timeline_until_horizon(mwpm, detection_events, horizon);
    ASSERT_EQ(status, TimelineStatus::TRUNCATED);

    auto exposed = surviving_exposed_regions(mwpm, detection_events);
    ASSERT_EQ(exposed.size(), 1u);
    ASSERT_EQ(blossom_nesting_depth(*exposed[0]), 1);

    // Exactly one member of the blossom has a dual equal to the horizon: that member is the base,
    // and no other member can stand in for it.
    std::vector<uint64_t> members;
    collect_defects_in_region(*exposed[0], mwpm.flooder.graph, members);
    std::sort(members.begin(), members.end());
    ASSERT_EQ(members, (std::vector<uint64_t>{0, 1, 2}));
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[2], horizon), horizon);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[0], horizon), 4);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[1], horizon), 2);

    auto harvest = harvest_to_obs(mwpm, detection_events);
    ASSERT_EQ(harvest.residual, (std::vector<uint64_t>{2}));
    ASSERT_EQ(harvest.residual_dual_sum, (std::vector<pm::total_weight_int>{horizon}));
    // The remaining even cycle is fully paired: r0 and r1, at their exact distance.
    ASSERT_EQ(harvest.committed_pairs_blossom_cycle, 1);
    ASSERT_EQ(harvest.committed.weight, 4);
    ASSERT_EQ(harvest.exposed_root_blossoms, 1);
    ASSERT_EQ(harvest.num_trees, 1);
    // y_2 + y_1 + y_0 + y_B = 5 + 1 + 3 + 1.
    ASSERT_EQ(harvest.dual_sum_at_truncation, 10);
    ASSERT_EQ(harvest.committed.weight + harvest.residual_dual_sum[0], harvest.dual_sum_at_truncation);
}

// Test 3, nested case: the base has to be recursed down through the blossom nesting. Returning any
// member of the outer blossom, or the base of the wrong sub-blossom, fails here.
TEST(TwoPhaseExposedBlossom, NestedBlossomExposesTheBaseOfItsBase) {
    auto mwpm = build_nested_blossom_fixture();
    std::vector<uint64_t> detection_events = {0, 1, 2, 3, 4};
    const pm::cumulative_time_int horizon = 13;

    auto status = process_timeline_until_horizon(mwpm, detection_events, horizon);
    ASSERT_EQ(status, TimelineStatus::TRUNCATED);

    auto exposed = surviving_exposed_regions(mwpm, detection_events);
    ASSERT_EQ(exposed.size(), 1u);
    ASSERT_EQ(blossom_nesting_depth(*exposed[0]), 2) << "the fixture did not produce a nested blossom";

    std::vector<uint64_t> members;
    collect_defects_in_region(*exposed[0], mwpm.flooder.graph, members);
    std::sort(members.begin(), members.end());
    ASSERT_EQ(members, (std::vector<uint64_t>{0, 1, 2, 3, 4}));
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[2], horizon), horizon);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[0], horizon), 11);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[1], horizon), 9);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[3], horizon), 3);
    ASSERT_EQ(nested_dual_sum(mwpm.flooder.graph.nodes[4], horizon), 5);

    auto harvest = harvest_to_obs(mwpm, detection_events);
    ASSERT_EQ(harvest.residual, (std::vector<uint64_t>{2}));
    ASSERT_EQ(harvest.residual_dual_sum, (std::vector<pm::total_weight_int>{horizon}));
    // One pair per level: (r0, r1) inside B, and (r3, r4) inside the outer blossom.
    ASSERT_EQ(harvest.committed_pairs_blossom_cycle, 2);
    ASSERT_EQ(harvest.committed.weight, 8);
    ASSERT_EQ(harvest.exposed_root_blossoms, 1);
    // y_0 + y_1 + y_2 + y_B + y_3 + y_4 + y_B2 = 3 + 1 + 5 + 6 + 1 + 3 + 2.
    ASSERT_EQ(harvest.dual_sum_at_truncation, 21);
}

// The match-edges flavour must expose the same defect and pair the same cycle.
TEST(TwoPhaseExposedBlossom, MatchEdgesFlavorAgreesOnTheNestedFixture) {
    auto mwpm = build_nested_blossom_fixture();
    std::vector<uint64_t> detection_events = {0, 1, 2, 3, 4};
    const pm::cumulative_time_int horizon = 13;

    ASSERT_EQ(process_timeline_until_horizon(mwpm, detection_events, horizon), TimelineStatus::TRUNCATED);
    std::vector<pm::CompressedEdge> match_edges;
    auto harvest = harvest_to_match_edges(mwpm, detection_events, match_edges);

    ASSERT_EQ(harvest.residual, (std::vector<uint64_t>{2}));
    std::vector<std::pair<int64_t, int64_t>> pairs;
    for (const auto& edge : match_edges) {
        int64_t from = edge.loc_from - mwpm.flooder.graph.nodes.data();
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)(edge.loc_to - mwpm.flooder.graph.nodes.data());
        pairs.push_back(from <= to ? std::make_pair(from, to) : std::make_pair(to, from));
    }
    std::sort(pairs.begin(), pairs.end());
    ASSERT_EQ(pairs, (std::vector<std::pair<int64_t, int64_t>>{{0, 1}, {3, 4}}));
}

// Test 3, corpus sweep. Real DEMs throw up exposed root blossoms of every shape; on all of them,
// exactly one member's dual reaches the horizon and that member is the one harvested. The failure
// mode this guards against is silent — a wrong member breaks the separation invariant downstream
// without anything crashing.
TEST(TwoPhaseExposedBlossom, CorpusExposedRootBlossomsAlwaysYieldTheirBase) {
    Corpus corpus = generate_surface_code_corpus(7, 7, 0.01, 250, 31);
    auto mwpm = corpus.to_mwpm();
    Harvester harvester;
    auto unit = median_edge_weight(mwpm.flooder.graph);

    size_t blossoms_seen = 0;
    size_t nested_blossoms_seen = 0;
    for (const auto& shot : corpus.shots) {
        for (int multiple : {1, 2, 3, 5, 8}) {
            pm::cumulative_time_int horizon = unit * multiple;
            auto status = process_timeline_until_horizon(mwpm, shot, horizon);
            if (status == TimelineStatus::COMPLETE) {
                harvester.harvest_to_obs(mwpm, shot);
                continue;
            }

            // Independently of the harvester, work out which defects the residual is allowed to
            // contain: per surviving tree, the members of its exposed region whose dual has
            // reached the horizon.
            //
            // Usually there is exactly one such member — the structural base. Ties do occur on a
            // degenerate graph like the surface code: a region that is matched and then pulled
            // back into a tree at the *same* instant spends zero time frozen, so its dual also
            // tracks the clock. Any member at the horizon is a legitimate exposed defect: the odd
            // cycle's edges are all tight, so removing any member leaves a tight perfect matching
            // on the rest, and the separation invariant only asks for `Y(u) == T`.
            std::vector<std::vector<uint64_t>> candidates_per_tree;
            for (pm::GraphFillRegion* exposed : surviving_exposed_regions(mwpm, shot)) {
                if (!exposed->blossom_children.empty()) {
                    blossoms_seen++;
                    if (blossom_nesting_depth(*exposed) >= 2)
                        nested_blossoms_seen++;
                }
                std::vector<uint64_t> members;
                collect_defects_in_region(*exposed, mwpm.flooder.graph, members);
                ASSERT_FALSE(members.empty());

                std::vector<uint64_t> candidates;
                for (uint64_t member : members) {
                    auto dual_sum = nested_dual_sum(mwpm.flooder.graph.nodes[member], horizon);
                    ASSERT_LE(dual_sum, horizon) << "a dual overshot the horizon";
                    if (dual_sum == horizon)
                        candidates.push_back(member);
                }
                ASSERT_FALSE(candidates.empty()) << "an exposed region had no member at the horizon";
                candidates_per_tree.push_back(std::move(candidates));
            }

            auto harvest = harvester.harvest_to_obs(mwpm, shot);
            ASSERT_EQ(harvest.residual.size(), candidates_per_tree.size());

            // Every harvested defect is a base of exactly one tree, and no tree is represented
            // twice: the residual really is one exposed defect per surviving tree.
            std::vector<bool> tree_claimed(candidates_per_tree.size(), false);
            for (uint64_t residual_defect : harvest.residual) {
                bool claimed = false;
                for (size_t tree = 0; tree < candidates_per_tree.size(); tree++) {
                    const auto& candidates = candidates_per_tree[tree];
                    if (tree_claimed[tree] ||
                        std::find(candidates.begin(), candidates.end(), residual_defect) == candidates.end())
                        continue;
                    tree_claimed[tree] = true;
                    claimed = true;
                    break;
                }
                ASSERT_TRUE(claimed) << "harvest returned " << residual_defect
                                     << ", which is not a base of any surviving tree, at T=" << horizon;
            }
        }
    }
    ASSERT_GT(blossoms_seen, 0u) << "no exposed root blossom was produced, so nothing was tested";
    ASSERT_GT(nested_blossoms_seen, 0u) << "no nested exposed root blossom was produced";
}
