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

/// §M2.9.7 — the harvest restructuring tests.
///
/// H1 is the deliverable: §M2.9's enumeration is bit-exact against M1.3's on the full corpus, both
/// flavours, with no tolerance. A failure here is a bug in M2.9, never a degeneracy — the two paths
/// read the *same* primal and dual state and apply the *same* commit policy, so there is no tie for
/// them to break differently.

#include <algorithm>
#include <map>

#include "gtest/gtest.h"

#include "pyrematching/sparse_blossom/flooder/graph_flooder.h"
#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

/// The committed pairing, in detector ids, canonicalised so that the *set* is compared rather than
/// the order the two enumerations happen to emit it in. `-1` is the boundary.
std::vector<std::pair<int64_t, int64_t>> canonical_pairs(
    const std::vector<pm::CompressedEdge>& match_edges, const pm::MatchingGraph& graph) {
    std::vector<std::pair<int64_t, int64_t>> pairs;
    pairs.reserve(match_edges.size());
    for (const pm::CompressedEdge& edge : match_edges) {
        int64_t from = (int64_t)(edge.loc_from - graph.nodes.data());
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)(edge.loc_to - graph.nodes.data());
        if (to >= 0 && to < from)
            std::swap(from, to);
        pairs.emplace_back(from, to);
    }
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::string describe(const std::vector<std::pair<int64_t, int64_t>>& pairs) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < pairs.size(); i++)
        out << (i ? ", " : "") << "(" << pairs[i].first << ", " << pairs[i].second << ")";
    out << "]";
    return out.str();
}

std::vector<Corpus> h1_corpora() {
    std::vector<Corpus> corpora;
    corpora.push_back(load_surface_code_d13(40));
    corpora.push_back(load_surface_code_d13_negative_weights(40));
    corpora.push_back(load_toric_code_d5(40));
    corpora.push_back(generate_surface_code_corpus(5, 5, 0.008, 40, 23));
    corpora.push_back(generate_surface_code_corpus(9, 9, 0.01, 40, 7));
    return corpora;
}

}  // namespace

// H1, obs flavour. Bit-exact against M1.3's landed harvest: residual vector, residual duals,
// committed obs bytes and weight, every pair counter, `num_trees`, `largest_tree_size`,
// `exposed_root_blossoms`, `dual_sum_at_truncation`. No tolerance.
TEST(TwoPhaseHarvestM29, H1ObsFlavourIsBitExactAgainstTheM1Oracle) {
    size_t truncated_shots = 0;
    for (const Corpus& corpus : h1_corpora()) {
        auto actual_mwpm = corpus.to_mwpm();
        auto oracle_mwpm = corpus.to_mwpm();
        Harvester actual;
        Harvester oracle;
        // `largest_tree_size` is a diagnostic on the new path and free on the old one, so the
        // comparison only means anything with diagnostics on.
        actual.collect_diagnostics = true;
        oracle.collect_diagnostics = true;
        oracle.use_legacy_enumeration = true;
        auto unit = median_edge_weight(actual_mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {0, 1, 2, 3, 5, 8}) {
                pm::cumulative_time_int horizon = unit * multiple;

                process_timeline_until_horizon(actual_mwpm, shot, horizon);
                HarvestResult got = actual.harvest_to_obs(actual_mwpm, shot);
                process_timeline_until_horizon(oracle_mwpm, shot, horizon);
                HarvestResult want = oracle.harvest_to_obs(oracle_mwpm, shot);

                ASSERT_TRUE(got.identical_output_to(want))
                    << corpus.name << " at T=" << horizon << ":" << got.describe_difference(want);
                if (!got.residual.empty())
                    truncated_shots++;
            }
        }
    }
    ASSERT_GT(truncated_shots, 0u) << "no shot ever truncated, so the tree path was never compared";
}

// H1, match-edges flavour. Same, plus the committed pair set itself, which only this flavour emits.
TEST(TwoPhaseHarvestM29, H1MatchEdgesFlavourIsBitExactAgainstTheM1Oracle) {
    size_t pairs_compared = 0;
    for (const Corpus& corpus : h1_corpora()) {
        auto actual_mwpm = corpus.to_mwpm();
        auto oracle_mwpm = corpus.to_mwpm();
        Harvester actual;
        Harvester oracle;
        actual.collect_diagnostics = true;
        oracle.collect_diagnostics = true;
        oracle.use_legacy_enumeration = true;
        auto unit = median_edge_weight(actual_mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {0, 1, 2, 3, 5, 8}) {
                pm::cumulative_time_int horizon = unit * multiple;

                std::vector<pm::CompressedEdge> got_edges;
                std::vector<pm::CompressedEdge> want_edges;
                process_timeline_until_horizon(actual_mwpm, shot, horizon);
                HarvestResult got = actual.harvest_to_match_edges(actual_mwpm, shot, got_edges);
                process_timeline_until_horizon(oracle_mwpm, shot, horizon);
                HarvestResult want = oracle.harvest_to_match_edges(oracle_mwpm, shot, want_edges);

                ASSERT_TRUE(got.identical_output_to(want))
                    << corpus.name << " at T=" << horizon << ":" << got.describe_difference(want);

                auto got_pairs = canonical_pairs(got_edges, actual_mwpm.flooder.graph);
                auto want_pairs = canonical_pairs(want_edges, oracle_mwpm.flooder.graph);
                ASSERT_EQ(got_pairs, want_pairs) << corpus.name << " at T=" << horizon << ": committed pair set moved, "
                                                 << describe(want_pairs) << " vs " << describe(got_pairs);
                pairs_compared += got_pairs.size();
            }
        }
    }
    ASSERT_GT(pairs_compared, 0u) << "no pair was ever committed, so nothing was compared";
}

// H2. The maintained enumeration equals the detection-event sweep's region set, on every shot.
//
// The debug build asserts this inline inside every harvest (§M2.9.2), which is the real guard; this
// test exists so that the guard is provably exercised across the corpus, and so that the equality
// is stated somewhere a release build can read.
TEST(TwoPhaseHarvestM29, H2LiveRegionListEqualsTheDetectionEventSweep) {
    size_t shots_checked = 0;
    for (const Corpus& corpus : h1_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 3, 6}) {
                process_timeline_until_horizon(mwpm, shot, unit * multiple);

                std::set<pm::GraphFillRegion*> from_live_list;
                for (pm::GraphFillRegion* region : mwpm.flooder.region_arena.live) {
                    if (region->blossom_parent == nullptr)
                        from_live_list.insert(region);
                }

                std::set<pm::GraphFillRegion*> from_sweep;
                const std::vector<uint64_t>* event_lists[2] = {
                    &shot, &mwpm.flooder.negative_weight_detection_events};
                for (const std::vector<uint64_t>* events : event_lists) {
                    for (uint64_t det : *events) {
                        if (det >= mwpm.flooder.graph.nodes.size())
                            continue;
                        const pm::DetectorNode& node = mwpm.flooder.graph.nodes[det];
                        if (node.region_that_arrived == nullptr)
                            continue;
                        from_sweep.insert(node.region_that_arrived_top);
                    }
                }

                ASSERT_EQ(from_live_list, from_sweep) << corpus.name << ": the maintained top-level region set is not "
                                                         "the one the detection-event sweep produces";
                shots_checked++;
                harvester.harvest_to_obs(mwpm, shot);
            }
        }
    }
    ASSERT_GT(shots_checked, 0u);
}

// H3. Committed and residual partition the shot's detection events, each classified exactly once,
// re-asserted against the new enumeration (debug invariant 3).
TEST(TwoPhaseHarvestM29, H3CommittedAndResidualStillPartitionTheDetectionEvents) {
    for (const Corpus& corpus : h1_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {0, 1, 3, 6}) {
                auto seeded = expected_seeded_detection_events(mwpm, shot);
                process_timeline_until_horizon(mwpm, shot, unit * multiple);
                std::vector<pm::CompressedEdge> match_edges;
                auto harvest = harvester.harvest_to_match_edges(mwpm, shot, match_edges);

                std::vector<uint64_t> classified = harvest.residual;
                for (const auto& edge : match_edges) {
                    classified.push_back((uint64_t)(edge.loc_from - mwpm.flooder.graph.nodes.data()));
                    if (edge.loc_to != nullptr)
                        classified.push_back((uint64_t)(edge.loc_to - mwpm.flooder.graph.nodes.data()));
                }
                std::sort(classified.begin(), classified.end());
                ASSERT_EQ(std::adjacent_find(classified.begin(), classified.end()), classified.end())
                    << corpus.name << ": a detection event was classified twice";
                ASSERT_EQ(classified, seeded) << corpus.name << ": committed + residual is not a partition";
                ASSERT_EQ(harvest.residual.size(), (size_t)harvest.num_trees);
            }
        }
    }
}

// H4. The separation invariant, unchanged from M1 test 2 and re-run against the new enumeration.
// Release blocker if it fires, and the base descent is the most likely cause.
TEST(TwoPhaseHarvestM29, H4SeparationInvariantSurvivesTheRestructuring) {
    size_t residual_defects_seen = 0;
    for (const Corpus& corpus : h1_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3, 5, 8}) {
                pm::cumulative_time_int horizon = unit * multiple;
                auto decoded = phase1_decode_to_obs(mwpm, shot, horizon, &harvester);
                for (size_t i = 0; i < decoded.harvest.residual.size(); i++) {
                    ASSERT_EQ(decoded.harvest.residual_dual_sum[i], horizon)
                        << corpus.name << ": Y(u) != T for residual defect " << decoded.harvest.residual[i];
                    residual_defects_seen++;
                }
            }
        }
    }
    ASSERT_GT(residual_defects_seen, 0u);
}

// H5. Base descent. M1.4's nested-blossom behaviour re-run against the mechanism §M2.9.4 selected
// (which is M1.4's recursion, unchanged): the extracted defect is the blossom's base, and it is the
// same defect the oracle path extracts.
TEST(TwoPhaseHarvestM29, H5ExposedRootBlossomsStillYieldTheirBase) {
    size_t blossom_roots_seen = 0;
    for (const Corpus& corpus : h1_corpora()) {
        auto actual_mwpm = corpus.to_mwpm();
        auto oracle_mwpm = corpus.to_mwpm();
        Harvester actual;
        Harvester oracle;
        actual.collect_diagnostics = true;
        oracle.collect_diagnostics = true;
        oracle.use_legacy_enumeration = true;
        auto unit = median_edge_weight(actual_mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3}) {
                pm::cumulative_time_int horizon = unit * multiple;

                process_timeline_until_horizon(actual_mwpm, shot, horizon);
                HarvestResult got = actual.harvest_to_obs(actual_mwpm, shot);
                process_timeline_until_horizon(oracle_mwpm, shot, horizon);
                HarvestResult want = oracle.harvest_to_obs(oracle_mwpm, shot);

                ASSERT_EQ(got.exposed_root_blossoms, want.exposed_root_blossoms) << corpus.name;
                ASSERT_EQ(got.residual, want.residual) << corpus.name << ": a root blossom exposed a different member";
                for (auto y : got.residual_dual_sum)
                    ASSERT_EQ(y, horizon) << corpus.name << ": the exposed member is not the base";
                blossom_roots_seen += (size_t)got.exposed_root_blossoms;
            }
        }
    }
    ASSERT_GT(blossom_roots_seen, 0u) << "no exposed root blossom ever occurred, so nothing was tested";
}

// H7. Reusability and determinism. A second shot decoded on a restructured-harvest instance matches
// a fresh instance exactly, and two runs of the same sequence are bit-identical.
TEST(TwoPhaseHarvestM29, H7ReusabilityAndDeterminism) {
    Corpus corpus = load_surface_code_d13(30);
    auto reused = corpus.to_mwpm();
    Harvester reused_harvester;
    reused_harvester.collect_diagnostics = true;
    auto unit = median_edge_weight(reused.flooder.graph);

    std::vector<HarvestResult> first_run;
    for (const auto& shot : corpus.shots) {
        process_timeline_until_horizon(reused, shot, unit * 2);
        first_run.push_back(reused_harvester.harvest_to_obs(reused, shot));
    }

    // Same sequence again on a second long-lived instance: bit-identical.
    auto reused_again = corpus.to_mwpm();
    Harvester again_harvester;
    again_harvester.collect_diagnostics = true;
    for (size_t i = 0; i < corpus.shots.size(); i++) {
        process_timeline_until_horizon(reused_again, corpus.shots[i], unit * 2);
        HarvestResult second = again_harvester.harvest_to_obs(reused_again, corpus.shots[i]);
        ASSERT_TRUE(second.identical_output_to(first_run[i]))
            << "shot " << i << " was not reproduced:" << second.describe_difference(first_run[i]);
    }

    // Every shot also matches a completely fresh instance, so nothing leaks between shots.
    for (size_t i = 0; i < corpus.shots.size(); i++) {
        auto fresh = corpus.to_mwpm();
        Harvester fresh_harvester;
        fresh_harvester.collect_diagnostics = true;
        process_timeline_until_horizon(fresh, corpus.shots[i], unit * 2);
        HarvestResult expected = fresh_harvester.harvest_to_obs(fresh, corpus.shots[i]);
        ASSERT_TRUE(first_run[i].identical_output_to(expected))
            << "shot " << i << " on a reused instance differs from a fresh one:"
            << first_run[i].describe_difference(expected);
    }
}

// H8. Fuzz: random horizons, random syndromes, small generated DEMs. H1 and H3 hold in all cases.
TEST(TwoPhaseHarvestM29, H8FuzzRandomHorizonsAndSyndromes) {
    std::mt19937_64 rng(20260803);
    size_t cases = 0;
    for (size_t distance : {3u, 5u, 7u}) {
        for (double noise : {0.005, 0.02, 0.05}) {
            Corpus corpus = generate_surface_code_corpus(distance, distance, noise, 25, rng());
            auto actual_mwpm = corpus.to_mwpm();
            auto oracle_mwpm = corpus.to_mwpm();
            Harvester actual;
            Harvester oracle;
            actual.collect_diagnostics = true;
            oracle.collect_diagnostics = true;
            oracle.use_legacy_enumeration = true;
            auto unit = median_edge_weight(actual_mwpm.flooder.graph);

            for (const auto& shot : corpus.shots) {
                pm::cumulative_time_int horizon = (pm::cumulative_time_int)(rng() % (6 * unit + 1));
                auto seeded = expected_seeded_detection_events(actual_mwpm, shot);

                std::vector<pm::CompressedEdge> got_edges;
                std::vector<pm::CompressedEdge> want_edges;
                process_timeline_until_horizon(actual_mwpm, shot, horizon);
                HarvestResult got = actual.harvest_to_match_edges(actual_mwpm, shot, got_edges);
                process_timeline_until_horizon(oracle_mwpm, shot, horizon);
                HarvestResult want = oracle.harvest_to_match_edges(oracle_mwpm, shot, want_edges);

                ASSERT_TRUE(got.identical_output_to(want))
                    << corpus.name << " at T=" << horizon << ":" << got.describe_difference(want);
                ASSERT_EQ(
                    canonical_pairs(got_edges, actual_mwpm.flooder.graph),
                    canonical_pairs(want_edges, oracle_mwpm.flooder.graph));

                std::vector<uint64_t> classified = got.residual;
                for (const auto& edge : got_edges) {
                    classified.push_back((uint64_t)(edge.loc_from - actual_mwpm.flooder.graph.nodes.data()));
                    if (edge.loc_to != nullptr)
                        classified.push_back((uint64_t)(edge.loc_to - actual_mwpm.flooder.graph.nodes.data()));
                }
                std::sort(classified.begin(), classified.end());
                ASSERT_EQ(classified, seeded);
                cases++;
            }
        }
    }
    ASSERT_GT(cases, 0u);
}

// §M2.9.6 measurement 1–3, exercised as a test so that the numbers the exit artifact reports are
// known to be filled in, and so that the diagnostics-off default really is free of them.
TEST(TwoPhaseHarvestM29, DiagnosticsAreFilledOnlyWhenAskedFor) {
    Corpus corpus = load_surface_code_d13(20);
    auto mwpm = corpus.to_mwpm();
    auto unit = median_edge_weight(mwpm.flooder.graph);

    Harvester quiet;
    Harvester loud;
    loud.collect_diagnostics = true;

    int loud_trees = 0;
    for (const auto& shot : corpus.shots) {
        process_timeline_until_horizon(mwpm, shot, unit);
        HarvestResult a = quiet.harvest_to_obs(mwpm, shot);
        ASSERT_EQ(a.largest_tree_size, 0) << "diagnostics leaked into the default path";
        ASSERT_EQ(a.max_blossom_nesting_depth, 0);
        ASSERT_EQ(a.harvest_dependent_depth, 0);

        process_timeline_until_horizon(mwpm, shot, unit);
        HarvestResult b = loud.harvest_to_obs(mwpm, shot);
        ASSERT_EQ(a.residual, b.residual) << "collect_diagnostics changed the output";
        ASSERT_EQ(a.committed, b.committed) << "collect_diagnostics changed the output";
        if (b.num_trees > 0) {
            ASSERT_GE(b.largest_tree_size, 1);
            ASSERT_GE(b.harvest_dependent_depth, 1);
            loud_trees++;
        }
    }
    ASSERT_GT(loud_trees, 0);
}
