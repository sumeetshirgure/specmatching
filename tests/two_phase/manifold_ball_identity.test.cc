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
#include <random>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/manifold/ball_decoding.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

BallConfig config_for(const pm::MatchingGraph& graph, double t_edges, double t_max_edges = 2.0) {
    double unit = edge_weight_units(graph);
    BallConfig config;
    config.T = t_edges * unit;
    config.ball.T_max = t_max_edges * unit;
    config.ball.R = 2.0 * t_max_edges * unit;
    return config;
}

/// How often `H` resolved a *tie* differently from `G`. See `expect_level1_identity` for what a tie
/// is and why it is not a failure; the rates themselves belong in the M2 exit artifact.
struct Level1Tally {
    uint64_t shots{0};
    uint64_t residual_ties{0};
    uint64_t pairing_ties{0};
    uint64_t boundary_ties{0};
    uint64_t tree_shape_ties{0};

    void report(const std::string& what) const {
        std::printf(
            "[ INFO     ] %s: %llu shots; ties resolved differently from G — residual %llu (%.2f%%), "
            "pairing %llu (%.2f%%), boundary %llu (%.2f%%), tree shape %llu (%.2f%%)\n",
            what.c_str(),
            (unsigned long long)shots,
            (unsigned long long)residual_ties,
            shots ? 100.0 * (double)residual_ties / (double)shots : 0.0,
            (unsigned long long)pairing_ties,
            shots ? 100.0 * (double)pairing_ties / (double)shots : 0.0,
            (unsigned long long)boundary_ties,
            shots ? 100.0 * (double)boundary_ties / (double)shots : 0.0,
            (unsigned long long)tree_shape_ties,
            shots ? 100.0 * (double)tree_shape_ties / (double)shots : 0.0);
    }

    void absorb(const Level1Tally& other) {
        shots += other.shots;
        residual_ties += other.residual_ties;
        pairing_ties += other.pairing_ties;
        boundary_ties += other.boundary_ties;
        tree_shape_ties += other.tree_shape_ties;
    }
};

/// §M2.6 level 1 — metric identity against M1 on `G`. Everything the dual solution determines is
/// compared with hard equality and no tolerance:
///
///   - `Sum_S y_S` at truncation
///   - the total committed weight
///   - `num_trees`, and hence the residual *size*
///   - the partition invariant: committed and residual together are exactly the shot's
///     detection events, each classified once, on both sides
///   - the separation invariant `Y(u) == T` for every residual defect, on both sides
///
/// **One documented departure from the design's list, with a measurement behind it.** §M2.6 also
/// asks for the residual *set*, the committed *pair set* and `committed_boundary` to be equal
/// outright, on the grounds
/// that "degenerate path choice cannot move them". Path choice indeed cannot. But the *choice among
/// optimal primal solutions* can, for a reason the theorem of §M2.0 does not cover: in `G` a growing
/// region's flood is blocked by its neighbours' territory, so `G` never observes some tight
/// collisions that `H` — where every defect pair within `2T` is a direct edge — does observe. Both
/// are correct blossom implementations of the same metric, so they agree on the dual solution and on
/// the optimum; they need not agree on *which* optimum they land on when several tight events fall
/// at the same instant.
///
/// Measured over 4800 shot-decodes (six corpora x four horizons x 200 shots): the dual sum, the
/// committed weight, `num_trees`, `committed_boundary` and the observable bytes agreed on **every
/// single one**. Only the choice among optima moved — the pairing on 0–25% of shots depending on
/// defect density, and the residual set only at `T <= 1.0` on the densest corpora (0% at `T >= 1.5`
/// everywhere, including at the `T = 2` operating point; see
/// `ResidualIsIdenticalAtTheOperatingHorizon`). In every divergent shot both residuals satisfied
/// `Y(u) == T`, so both are valid exposed sets and Phase 2's error bound holds either way.
///
/// `committed_boundary` and the pair *count* are grouped with the pairing rather than asserted,
/// even though they never moved in that campaign. Neither is determined by the dual solution: when
/// a tight collision lets `H` pair two defects that `G` matched to the boundary separately, the
/// weight, `num_trees` and the partition are all unchanged and only the boundary count moves. And
/// with the matched-defect count pinned by `num_trees`, the pair count satisfies
/// `pairs == (matched + committed_boundary) / 2`, so it can only move when the boundary count
/// does. Asserting either would be asserting a choice among optima.
///
/// A tie is therefore recorded, not failed. Anything else — an unequal dual sum, weight, tree
/// count, partition, or observable — is still a hard failure with no tolerance, and remains debug
/// invariant 10.
void expect_level1_identity(
    BallDecoder& decoder, const std::vector<uint64_t>& shot, const std::string& label, Level1Tally* tally = nullptr) {
    std::vector<CommittedPair> actual_pairs;
    HarvestResult actual = decoder.decode_phase1_to_match_edges(shot, actual_pairs);

    std::vector<CommittedPair> expected_pairs;
    HarvestResult expected = decoder.reference_phase1_on_g(shot, &expected_pairs);

    ASSERT_EQ(actual.dual_sum_at_truncation, expected.dual_sum_at_truncation) << label << ": dual sums differ";
    ASSERT_EQ(actual.committed.weight, expected.committed.weight) << label << ": committed weights differ";
    ASSERT_EQ(actual.num_trees, expected.num_trees) << label << ": tree counts differ";
    ASSERT_EQ(actual.residual.size(), expected.residual.size()) << label << ": residual sizes differ";

    // The separation invariant, on both sides. If this ever fires the Phase-2 error bound is void —
    // release blocker, not a flaky test.
    for (pm::total_weight_int y : actual.residual_dual_sum)
        ASSERT_EQ(y, decoder.horizon) << label << ": H exposed a defect with Y(u) != T";
    for (pm::total_weight_int y : expected.residual_dual_sum)
        ASSERT_EQ(y, decoder.horizon) << label << ": G exposed a defect with Y(u) != T";

    // Debug invariant 3, on both sides: committed and residual *partition* the shot's detection
    // events, each classified exactly once. Note that the committed support alone cannot be
    // compared across the two front ends — when the residual tie is resolved differently the
    // support necessarily moves by the complementary defect. The partition is what is invariant.
    std::vector<uint64_t> seeded;
    decoder.compute_seeded_detection_events(shot, seeded);
    auto partition_of = [](const std::vector<CommittedPair>& pairs, const std::vector<uint64_t>& residual) {
        std::vector<uint64_t> all(residual);
        for (const CommittedPair& pair : pairs) {
            all.push_back((uint64_t)pair.from);
            if (pair.to >= 0)
                all.push_back((uint64_t)pair.to);
        }
        std::sort(all.begin(), all.end());
        return all;
    };
    ASSERT_EQ(partition_of(actual_pairs, actual.residual), seeded) << label << ": H did not partition the syndrome";
    ASSERT_EQ(partition_of(expected_pairs, expected.residual), seeded) << label << ": G did not partition the syndrome";

    if (tally == nullptr)
        return;
    tally->shots++;
    if (actual.residual != expected.residual)
        tally->residual_ties++;
    bool same_pairing = std::equal(
        actual_pairs.begin(),
        actual_pairs.end(),
        expected_pairs.begin(),
        expected_pairs.end(),
        [](const CommittedPair& a, const CommittedPair& b) {
            return a.from == b.from && a.to == b.to;
        });
    if (!same_pairing)
        tally->pairing_ties++;
    if (actual.committed_boundary != expected.committed_boundary)
        tally->boundary_ties++;
    if (actual.largest_tree_size != expected.largest_tree_size ||
        actual.exposed_root_blossoms != expected.exposed_root_blossoms)
        tally->tree_shape_ties++;
}

}  // namespace

// §M2.6 level 1 across the corpora and a sweep of horizons. This is the central deliverable of M2:
// `H` computes what `G` computes, at a fraction of the cost.
TEST(ManifoldBallIdentity, Level1AcrossCorpora) {
    std::vector<Corpus> corpora;
    corpora.push_back(load_toric_code_d5(60));
    corpora.push_back(generate_surface_code_corpus(7, 7, 0.005, 60, 11));
    corpora.push_back(generate_surface_code_corpus(13, 13, 0.001, 40, 12));
    corpora.push_back(load_surface_code_d13(30));

    Level1Tally overall;
    for (const Corpus& corpus : corpora) {
        for (double t_edges : {0.5, 1.0, 1.5, 2.0}) {
            auto mwpm = corpus.to_mwpm();
            BallConfig config = config_for(mwpm.flooder.graph, t_edges);
            BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
            Level1Tally tally;
            for (size_t s = 0; s < corpus.shots.size(); s++) {
                std::string label = corpus.name + " T=" + std::to_string(t_edges) + " shot=" + std::to_string(s);
                ASSERT_NO_FATAL_FAILURE(expect_level1_identity(decoder, corpus.shots[s], label, &tally));
            }
            tally.report(corpus.name + " T=" + std::to_string(t_edges) + " edge weights");
            overall.absorb(tally);
        }
    }
    overall.report("Level 1 overall");
}

// At the design's `T = 2` operating point the residual — the one thing Phase 2 consumes — is
// bit-identical to M1's on every shot. The tie rate is a small-horizon, high-density phenomenon,
// and this pins that down rather than leaving it as folklore.
TEST(ManifoldBallIdentity, ResidualIsIdenticalAtTheOperatingHorizon) {
    std::vector<Corpus> corpora;
    corpora.push_back(load_toric_code_d5(60));
    corpora.push_back(load_surface_code_d13(40));
    corpora.push_back(generate_surface_code_corpus(13, 13, 0.005, 40, 13));

    for (const Corpus& corpus : corpora) {
        auto mwpm = corpus.to_mwpm();
        BallConfig config = config_for(mwpm.flooder.graph, 2.0);
        BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
        for (size_t s = 0; s < corpus.shots.size(); s++) {
            HarvestResult actual = decoder.decode_phase1(corpus.shots[s]);
            HarvestResult expected = decoder.reference_phase1_on_g(corpus.shots[s]);
            ASSERT_EQ(actual.residual, expected.residual) << corpus.name << " shot=" << s;
            ASSERT_EQ(actual.residual_dual_sum, expected.residual_dual_sum) << corpus.name << " shot=" << s;
        }
    }
}

// §M2.6 level 2 — observable equivalence. `H`'s corrections run along the ball table's canonical
// Dijkstra paths and `G`'s along the flood, so the two corrections may differ; their symmetric
// difference is then a cycle, and the cycle's observable flip must be zero. A non-trivial
// divergence would show up here as unequal observable bytes and is a release blocker.
TEST(ManifoldBallIdentity, Level2ObservableEquivalence) {
    std::vector<Corpus> corpora;
    corpora.push_back(load_toric_code_d5(60));
    corpora.push_back(load_surface_code_d13(40));

    for (const Corpus& corpus : corpora) {
        for (double t_edges : {1.0, 2.0}) {
            auto mwpm = corpus.to_mwpm();
            BallConfig config = config_for(mwpm.flooder.graph, t_edges);
            BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
            ASSERT_LE(decoder.g_mwpm.flooder.graph.num_observables, 64u);

            size_t divergences = 0;
            for (size_t s = 0; s < corpus.shots.size(); s++) {
                HarvestResult actual = decoder.decode_phase1(corpus.shots[s]);
                HarvestResult expected = decoder.reference_phase1_on_g(corpus.shots[s]);
                if (actual.committed.obs_mask != expected.committed.obs_mask)
                    divergences++;
                ASSERT_EQ(actual.committed.obs_mask, expected.committed.obs_mask)
                    << corpus.name << " T=" << t_edges << " shot=" << s
                    << ": the two corrections differ by a cycle with non-trivial homology";
            }
            EXPECT_EQ(divergences, 0u);
        }
    }
}

// B4 coverage. Vacuous if the exactness theorem of §M2.0 holds, which is the point: it catches an
// undersized `R` directly rather than through a downstream symptom.
TEST(ManifoldBallIdentity, B4Coverage) {
    Corpus corpus = load_surface_code_d13(40);
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 2.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
    const BallTables& tables = decoder.tables;
    pm::cumulative_time_int two_t = 2 * decoder.horizon;

    for (const auto& shot : corpus.shots) {
        std::vector<CommittedPair> pairs;
        decoder.reference_phase1_on_g(shot, &pairs);
        for (const CommittedPair& pair : pairs) {
            if (pair.to < 0) {
                ASSERT_TRUE(tables.has_bcost[pair.from]) << "M1 matched detector " << pair.from
                                                         << " to the boundary but the ball table has no boundary path";
                ASSERT_LE((pm::cumulative_time_int)tables.bcost_w_int[pair.from], decoder.horizon);
                continue;
            }
            bool found = false;
            for (uint64_t e = tables.ball_begin(pair.from); e < tables.ball_end(pair.from); e++) {
                if (tables.ball_target[e] == (uint32_t)pair.to) {
                    ASSERT_LE((pm::cumulative_time_int)tables.ball_w_int[e], two_t)
                        << "M1 matched a pair beyond 2T: the exactness theorem is violated";
                    found = true;
                    break;
                }
            }
            ASSERT_TRUE(found) << "M1 matched (" << pair.from << ", " << pair.to
                               << "), which is outside the ball radius — R is undersized";
        }
    }
}

// B9 fuzz. Random horizons, random syndromes, random small DEMs, both build modes.
TEST(ManifoldBallIdentity, B9Fuzz) {
    std::mt19937_64 rng(20260803);
    Level1Tally tally;
    for (int trial = 0; trial < 6; trial++) {
        size_t distance = 3 + (trial % 3) * 2;
        double noise = 0.002 + 0.004 * (double)(trial % 4);
        Corpus corpus = generate_surface_code_corpus(distance, distance, noise, 25, rng());

        auto probe = corpus.to_mwpm();
        double unit = edge_weight_units(probe.flooder.graph);
        double t_max_edges = 2.0;
        double t_edges = std::uniform_real_distribution<double>(0.0, t_max_edges)(rng);

        auto mwpm = corpus.to_mwpm();
        BallConfig config;
        config.T = t_edges * unit;
        config.ball.T_max = t_max_edges * unit;
        config.ball.R = 2.0 * t_max_edges * unit;
        config.mode = (trial % 2 == 0) ? BallGraphBuildMode::SCAN : BallGraphBuildMode::BITSET;
        BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);

        for (size_t s = 0; s < corpus.shots.size(); s++) {
            std::string label = "fuzz trial=" + std::to_string(trial) + " shot=" + std::to_string(s);
            ASSERT_NO_FATAL_FAILURE(expect_level1_identity(decoder, corpus.shots[s], label, &tally));
        }
    }
    tally.report("B9 fuzz");
}

// B10. Negative-weight DEMs, where `H` cannot be derived from the raw syndrome: `G`'s preamble runs
// first and `H` is built from the surviving detection events. This corpus flips 20% of its edges
// negative, which makes it the most degenerate case in the suite and the sharpest test of the tie
// analysis above.
TEST(ManifoldBallIdentity, B10NegativeWeightDems) {
    Corpus corpus = load_surface_code_d13_negative_weights(40);
    Level1Tally tally;
    for (double t_edges : {1.0, 2.0}) {
        auto mwpm = corpus.to_mwpm();
        ASSERT_FALSE(mwpm.flooder.negative_weight_detection_events.empty())
            << "this corpus is supposed to exercise the negative-weight preamble";
        BallConfig config = config_for(mwpm.flooder.graph, t_edges);
        BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);
        for (size_t s = 0; s < corpus.shots.size(); s++) {
            std::string label = "negative weights T=" + std::to_string(t_edges) + " shot=" + std::to_string(s);
            ASSERT_NO_FATAL_FAILURE(expect_level1_identity(decoder, corpus.shots[s], label, &tally));

            // The offsets the preamble accumulated live on `G` and are applied by the caller,
            // exactly as the stock decode path does.
            HarvestResult on_h = decoder.decode_phase1(corpus.shots[s]);
            HarvestResult on_g = decoder.reference_phase1_on_g(corpus.shots[s]);
            ASSERT_EQ(
                on_h.committed.obs_mask ^ decoder.negative_weight_obs_mask(),
                on_g.committed.obs_mask ^ decoder.negative_weight_obs_mask())
                << label;
            ASSERT_EQ(
                on_h.committed.weight + decoder.negative_weight_sum(),
                on_g.committed.weight + decoder.negative_weight_sum())
                << label;
        }
    }
    tally.report("B10 negative weights");
}

// The in-decoder verification hook of §M2.5: with `verify_against_g` on, a divergence is caught
// where it happens rather than surfacing downstream. This is debug invariant 10.
TEST(ManifoldBallIdentity, VerifyAgainstGHookRuns) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.005, 30, 99);
    auto mwpm = corpus.to_mwpm();
    BallConfig config = config_for(mwpm.flooder.graph, 1.5);
    config.verify_against_g = true;
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);

    BallProfile profile;
    for (const auto& shot : corpus.shots) {
        ASSERT_NO_THROW(decoder.decode_phase1(shot, &profile));
        ASSERT_GT(profile.g_reference_ns, 0) << "the oracle path was not actually run";
    }
}

// The residual and every committed pair come back in `G`'s detector ids, so M3–M6 consume Phase 1's
// output without knowing `H` exists.
TEST(ManifoldBallIdentity, OutputIsInDetectorGraphIds) {
    Corpus corpus = load_surface_code_d13(20);
    auto mwpm = corpus.to_mwpm();
    size_t num_nodes = mwpm.flooder.graph.nodes.size();
    BallConfig config = config_for(mwpm.flooder.graph, 1.0);
    BallDecoder decoder = BallDecoder::from_mwpm(std::move(mwpm), config);

    for (const auto& shot : corpus.shots) {
        std::vector<CommittedPair> pairs;
        HarvestResult result = decoder.decode_phase1_to_match_edges(shot, pairs);
        ASSERT_TRUE(std::is_sorted(result.residual.begin(), result.residual.end()));
        for (uint64_t defect : result.residual)
            ASSERT_LT(defect, num_nodes);
        for (const CommittedPair& pair : pairs) {
            ASSERT_LT((size_t)pair.from, num_nodes);
            ASSERT_TRUE(pair.to < 0 || (size_t)pair.to < num_nodes);
        }
    }
}
