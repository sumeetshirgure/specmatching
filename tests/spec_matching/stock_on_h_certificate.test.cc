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

/// §M7.6 — the stock-on-`H` front end and its max-dual certificate.
///
/// The oracle here is **stock exact decode on `G`**, not M1's truncated harvest: this path produces
/// no truncated intermediate state to reproduce, and a certified shot is a global MWPM by the
/// theorem of §M7.0. So every level is an equality against the inherited decoder, and a failure is a
/// certificate bug rather than an accuracy tradeoff. Do not add a tolerance to any of it.
///
/// C3 (forced escalation) is written first and deliberately: it is the cheapest end-to-end proof
/// that the plumbing is right, and it holds whatever the certificate does.

#include <cstdio>
#include <map>
#include <random>

#include "gtest/gtest.h"

#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/spec_matching/certificate/max_dual.h"
#include "pyrematching/spec_matching/driver/spec_matching_decoding.h"
#include "tests/spec_matching/spec_matching_test_util.h"

using namespace pm::spec_matching;
using namespace pm::spec_matching::test;

namespace {

SpecMatchingConfig stock_config_for(const pm::MatchingGraph& graph, double t_edges, double t_max_edges = 2.0) {
    double unit = edge_weight_units(graph);
    SpecMatchingConfig config;
    config.T = t_edges * unit;
    config.ball.T_max = t_max_edges * unit;
    config.ball.R = 2.0 * t_max_edges * unit;
    config.stock_on_h = true;
    return config;
}

struct Answer {
    std::vector<uint8_t> obs;
    pm::total_weight_int weight{0};
};

Answer stock_decode(pm::Mwpm& mwpm, const std::vector<uint64_t>& dets, size_t num_observables) {
    Answer answer;
    answer.obs.assign(num_observables, 0);
    pm::decode_detection_events(mwpm, dets, answer.obs.data(), answer.weight, false);
    return answer;
}

Answer spec_matching_decode(SpecMatchingDecoder& decoder, const std::vector<uint64_t>& dets, SpecMatchingProfile* prof = nullptr) {
    Answer answer;
    answer.obs.assign(decoder.num_observables, 0);
    decoder.decode_to_obs(dets, answer.obs.data(), answer.weight, prof);
    return answer;
}

using Pair = std::pair<int64_t, int64_t>;

Pair canonical(int64_t a, int64_t b) {
    // The boundary end always sits second, so a boundary match has one canonical form.
    if (b >= 0 && b < a)
        std::swap(a, b);
    return {a, b};
}

/// Stock exact decode's **matched pairs** on `G`, as sorted `(min, max)` detector-id pairs with
/// `-1` for the boundary. This is the level-1 oracle.
///
/// Written out here rather than routed through `decode_detection_events_to_match_edges`, which
/// refuses negative-weight graphs outright — C8 needs exactly those. The extraction mirrors
/// `decode_detection_events_for_up_to_64_observables`: both the shot's own detection events and the
/// graph's negative-weight ones, with `region_that_arrived` re-read each time so that a pair already
/// consumed from its other endpoint is not extracted twice.
std::vector<Pair> stock_pairs_on_g(pm::Mwpm& mwpm, const std::vector<uint64_t>& dets) {
    EXPECT_EQ(pm::process_timeline_until_completion_or_report(mwpm, dets), pm::CompletionStatus::COMPLETE);
    auto& nodes = mwpm.flooder.graph.nodes;
    mwpm.flooder.match_edges.clear();
    const std::vector<uint64_t>* event_lists[2] = {&dets, &mwpm.flooder.negative_weight_detection_events};
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det < nodes.size() && nodes[det].region_that_arrived != nullptr)
                mwpm.shatter_blossom_and_extract_match_edges(
                    nodes[det].region_that_arrived_top, mwpm.flooder.match_edges);
        }
    }

    const pm::DetectorNode* base = nodes.data();
    std::vector<Pair> pairs;
    pairs.reserve(mwpm.flooder.match_edges.size());
    for (const pm::CompressedEdge& edge : mwpm.flooder.match_edges) {
        int64_t from = (int64_t)(edge.loc_from - base);
        int64_t to = edge.loc_to == nullptr ? -1 : (int64_t)(edge.loc_to - base);
        pairs.push_back(canonical(from, to));
    }
    std::sort(pairs.begin(), pairs.end());
    mwpm.flooder.match_edges.clear();
    return pairs;
}

std::vector<Pair> to_pairs(const std::vector<CommittedPair>& committed) {
    std::vector<Pair> pairs;
    pairs.reserve(committed.size());
    for (const CommittedPair& pair : committed)
        pairs.push_back(canonical(pair.from, pair.to));
    std::sort(pairs.begin(), pairs.end());
    return pairs;
}

std::string describe(const std::vector<Pair>& pairs) {
    std::string text = "[";
    for (size_t i = 0; i < pairs.size(); i++) {
        text += (i ? ", (" : "(") + std::to_string(pairs[i].first) + ", " + std::to_string(pairs[i].second) + ")";
    }
    return text + "]";
}

// ---------------------------------------------------------------------------------------------
// The hand-built fixture behind C1 and C2.
//
// A star with centre `x = 4` and one shortcut, all weights in DEM float units:
//
//     a=0 --1-- x --1-- b=1          d(a, b)  =  2
//     c=2 --7-- x --7-- d=3          d(a, c)  =  d(a, d) = d(b, c) = d(b, d) = 8
//     c=2 -------10------ d=3        d(c, d)  = 10  (the shortcut beats 14 through x)
//
// At `T = 4` the ball filter admits every pair within `2T = 8` and excludes `(c, d)` at 10. So on
// the syndrome `{a, b, c, d}`:
//
//  - the true optimum is `(a,b) + (c,d) = 2 + 10 = 12`, and it needs the edge `H` does not have;
//  - `H` nonetheless *has* a perfect matching — `(a,c) + (b,d) = 16` — so the run **completes** and
//    nothing throws. This is §M7.0's "catching the throw is not sufficient" case, in the flesh.
//  - the certificate catches it: `Y(a) + Y(b) <= d(a,b) = 2` forces the duals onto the far pair, so
//    `Y(c), Y(d) >= 7 > T = 4` and the shot escalates.
//
// On the syndrome `{c, d}` the same graph gives C2: the only partner is at 10 > 2T and the graph
// has no boundary, so `H` has two nodes and no edges at all.
constexpr size_t FIXTURE_NODES = 5;
constexpr size_t FIXTURE_OBSERVABLES = 5;
constexpr double FIXTURE_T = 4.0;
constexpr double FIXTURE_OPTIMAL_WEIGHT = 12.0;
constexpr double FIXTURE_H_MATCHING_WEIGHT = 16.0;

pm::UserGraph fixture_graph() {
    pm::UserGraph graph(FIXTURE_NODES, FIXTURE_OBSERVABLES);
    graph.add_or_merge_edge(0, 4, {0}, 1.0, -1);
    graph.add_or_merge_edge(1, 4, {1}, 1.0, -1);
    graph.add_or_merge_edge(2, 4, {2}, 7.0, -1);
    graph.add_or_merge_edge(3, 4, {3}, 7.0, -1);
    graph.add_or_merge_edge(2, 3, {4}, 10.0, -1);
    return graph;
}

BallConfig fixture_config() {
    BallConfig config;
    config.T = FIXTURE_T;
    config.ball.T_max = FIXTURE_T;
    config.ball.R = 2 * FIXTURE_T;
    config.stock_on_h = true;
    return config;
}

/// The fixture's weights in the discretised metric the decoder reports in.
pm::total_weight_int in_time_units(const pm::Mwpm& mwpm, double weight) {
    return (pm::total_weight_int)to_time_units(weight, mwpm.flooder.graph.normalising_constant);
}

Corpus small_corpus(size_t shots = 60) {
    return generate_surface_code_corpus(5, 5, 0.01, shots, 20260901);
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// §M7.10 step 1 — the vendored entry point itself, in isolation from everything built on it.
// `process_timeline_until_completion_or_report` must return `NO_PERFECT_MATCHING` exactly where
// stock throws, and `COMPLETE` otherwise, and it must leave the teardown to its caller.
TEST(StockOnHCertificate, CompletionOrReportMatchesWhereStockThrows) {
    // The fixture graph has no boundary at all, so an odd-parity syndrome on it admits no perfect
    // matching — which is the case stock's own entry point raises `std::invalid_argument` for.
    pm::UserGraph graph = fixture_graph();
    pm::Mwpm reporting = graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false);
    pm::Mwpm throwing = graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false);

    std::vector<uint64_t> odd = {0};
    ASSERT_EQ(
        pm::process_timeline_until_completion_or_report(reporting, odd), pm::CompletionStatus::NO_PERFECT_MATCHING);
    // Reported, not thrown, and the instance is left standing — the caller owns the teardown,
    // because M7's certificate has to read the terminal state before anything is torn down.
    ASSERT_FALSE(reporting.flooder.region_arena.live.empty());
    reporting.reset();
    {
        pm::total_weight_int weight = 0;
        std::vector<uint8_t> obs(throwing.flooder.graph.num_observables, 0);
        ASSERT_THROW(pm::decode_detection_events(throwing, odd, obs.data(), weight, false), std::invalid_argument);
    }

    // And `COMPLETE` wherever stock does not throw, over a real corpus.
    Corpus corpus = small_corpus(60);
    pm::Mwpm mwpm = corpus.to_mwpm();
    for (const auto& shot : corpus.shots) {
        ASSERT_EQ(pm::process_timeline_until_completion_or_report(mwpm, shot), pm::CompletionStatus::COMPLETE);
        mwpm.reset();
    }
}

// ---------------------------------------------------------------------------------------------
// C3 — forced escalation. At a horizon small enough that nothing certifies, the decoder is stock,
// bit for bit, on the whole corpus. Written first: it exercises every stage of the new front end
// (ball filter, stock solve, certificate, escalation, teardown, reuse) and depends on none of the
// certificate's cleverness being right, only on its plumbing.
TEST(StockOnHCertificate, C3ForcedEscalationIsBitIdenticalToStock) {
    for (Corpus corpus : {load_surface_code_d13(40), load_surface_code_d13_negative_weights(40)}) {
        auto reference_mwpm = corpus.to_mwpm();
        // A quarter of one lattice edge: no pair and no boundary survives the `2T` filter, so every
        // shot with a defect in it escalates.
        SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, 0.25, 0.25);
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

        size_t escalated = 0;
        size_t with_defects = 0;
        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << corpus.name;
            ASSERT_EQ(actual.weight, expected.weight) << corpus.name;
            escalated += profile.escalated ? 1 : 0;
            if (!decoder.ball->arena.graph.h_to_det.empty())
                with_defects++;
            if (profile.escalated)
                ASSERT_EQ(profile.certified, 0);
        }
        ASSERT_GT(escalated, 0u) << corpus.name << ": the horizon was not small enough to force escalation";
        ASSERT_EQ(escalated, with_defects) << corpus.name << ": a shot with defects certified at T = 0.25 edges";
    }
}

// ---------------------------------------------------------------------------------------------
// Levels 1 and 2 — the identity the whole front end rests on. Matched pairs, total weight and
// observable bytes against stock exact decode on `G`, on **every** shot, certified or escalated, at
// several horizons. No tolerance.
//
// Level 2's licence to differ is narrower than it looks: `H` derives a correction from the ball
// table's canonical path where `G`'s flood may pick another minimum-weight path, so the two
// *corrections* can differ — but only by a cycle, and only a homologically trivial one, whose
// observable flip is zero. So the obs bytes still have to agree, and a divergence is a release
// blocker rather than a tolerated tie.
TEST(StockOnHCertificate, Level1And2AgainstStockOnG) {
    for (Corpus corpus :
         {load_surface_code_d13(40), load_surface_code_d13_negative_weights(40), load_toric_code_d5(40)}) {
        auto reference_mwpm = corpus.to_mwpm();
        auto pair_oracle_mwpm = corpus.to_mwpm();
        for (double t_edges : {0.75, 1.5, 2.0}) {
            SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, t_edges);
            auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

            BallConfig ball_config;
            ball_config.T = config.T;
            ball_config.ball = config.ball;
            ball_config.stock_on_h = true;
            auto pairs_decoder = BallDecoder::from_detector_error_model(corpus.dem, ball_config, NUM_DISTINCT_WEIGHTS);

            size_t certified = 0;
            size_t mask_divergences = 0;
            size_t pairing_ties = 0;
            std::vector<CommittedPair> committed;
            for (const auto& shot : corpus.shots) {
                SpecMatchingProfile profile;
                Answer actual = spec_matching_decode(decoder, shot, &profile);
                Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);

                // Level 1, weight half: hard equality on every shot.
                ASSERT_EQ(actual.weight, expected.weight)
                    << corpus.name << " at T = " << t_edges << " edges: a certified shot's weight is an exact MWPM "
                    << "weight, and an escalated shot's is stock's own";
                // Level 2.
                if (actual.obs != expected.obs)
                    mask_divergences++;
                ASSERT_EQ(actual.obs, expected.obs)
                    << corpus.name << " at T = " << t_edges
                    << " edges: the corrections differ by a cycle that flips an observable, which is a "
                       "non-trivial homology and a release blocker";

                // Level 1, pair half. The pairs come off a second decoder driven at the ball layer,
                // because the obs flavour never materialises them; both decoders are deterministic
                // and identically configured, so they take the same branch on the same shot.
                Phase1Outcome outcome = pairs_decoder.decode_phase1_production_to_match_edges(shot, committed);
                ASSERT_EQ(outcome.status == TimelineStatus::TRUNCATED, profile.escalated)
                    << "the two decoders disagreed on whether the shot escalates";
                if (outcome.status == TimelineStatus::COMPLETE) {
                    certified++;
                    std::vector<Pair> actual_pairs = to_pairs(committed);
                    std::vector<Pair> expected_pairs = stock_pairs_on_g(pair_oracle_mwpm, shot);
                    if (actual_pairs != expected_pairs) {
                        // A **tie**, not a divergence, and §M7.6's level 1 does not anticipate it —
                        // but M2 already measured it and §M2.6's own level 1 records exactly this.
                        // When the optimum is degenerate there is more than one minimum-weight
                        // perfect matching, and `H` and `G` need not land on the same one: in `G` a
                        // growing region's flood is blocked by its neighbours' territory, so `G`
                        // never observes some tight collisions that `H` — where every pair within
                        // `2T` is a direct edge — does. Both are correct blossom implementations of
                        // the same metric.
                        //
                        // What the certificate actually claims is that `H`'s matching is *a* global
                        // MWPM, and that claim is checked above and not weakened here: the weight
                        // assertion is hard, `H`'s weight is the sum of the ball tables' exact
                        // `d_G` values over its own pairs, and it equals stock's optimum. So a
                        // differing pairing at equal weight is a second optimum, and a differing
                        // pairing at a different weight would already have failed.
                        //
                        // What is still asserted here is the part a tie cannot move: the pairing
                        // covers the same defects, exactly once each, and has the same size.
                        pairing_ties++;
                        ASSERT_EQ(actual_pairs.size(), expected_pairs.size())
                            << corpus.name << " at T = " << t_edges << " edges: " << describe(actual_pairs) << " vs "
                            << describe(expected_pairs);
                        std::vector<int64_t> actual_support;
                        std::vector<int64_t> expected_support;
                        for (const Pair& pair : actual_pairs) {
                            actual_support.push_back(pair.first);
                            if (pair.second >= 0)
                                actual_support.push_back(pair.second);
                        }
                        for (const Pair& pair : expected_pairs) {
                            expected_support.push_back(pair.first);
                            if (pair.second >= 0)
                                expected_support.push_back(pair.second);
                        }
                        std::sort(actual_support.begin(), actual_support.end());
                        std::sort(expected_support.begin(), expected_support.end());
                        ASSERT_EQ(actual_support, expected_support)
                            << corpus.name << " at T = " << t_edges
                            << " edges: the two pairings do not even cover the same defects, so this is not a tie";
                    }
                }
            }
            std::printf(
                "[ INFO     ] %s at T = %.2f edges: %zu/%zu certified, %zu mask divergences, %zu pairing ties\n",
                corpus.name.c_str(),
                t_edges,
                certified,
                corpus.shots.size(),
                mask_divergences,
                pairing_ties);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// C1 — certificate soundness, on the fixture that makes the point: `H` completes with a
// **suboptimal** matching, so nothing throws, and only the dual test catches it. This is the case
// invariant 4 guards, and the reason the scheme is a certificate rather than an exception handler.
TEST(StockOnHCertificate, C1CompletingButOverHorizonEscalates) {
    pm::UserGraph graph = fixture_graph();
    pm::Mwpm reference = graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false);
    auto decoder = BallDecoder::from_mwpm(graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false), fixture_config());

    std::vector<uint64_t> shot = {0, 1, 2, 3};
    std::vector<CommittedPair> committed;
    BallProfile profile;
    Phase1Outcome outcome = decoder.decode_phase1_production_to_match_edges(shot, committed, &profile);

    // It completed — that is the whole point of the fixture — and it was rejected by the dual, not
    // by an exception.
    ASSERT_EQ(profile.h_no_perfect_matching, 0) << "the fixture is supposed to admit a perfect matching in H";
    ASSERT_EQ(profile.certified, 0);
    ASSERT_EQ(outcome.status, TimelineStatus::TRUNCATED) << "invariant 4: a completing-but-over-T H-matching escaped";
    ASSERT_GT(profile.max_dual_at_completion, decoder.horizon)
        << "the certificate rejected the shot for the wrong reason";

    // And escalation is not cosmetic: stock's answer is strictly better than the one `H` completed
    // with, so returning `H`'s would have been a silent wrong answer of exactly the size the
    // fixture was built to expose.
    Answer expected = stock_decode(reference, shot, reference.flooder.graph.num_observables);
    ASSERT_EQ(expected.weight, in_time_units(reference, FIXTURE_OPTIMAL_WEIGHT));
    ASSERT_LT(expected.weight, in_time_units(reference, FIXTURE_H_MATCHING_WEIGHT));
}

// The same soundness statement swept over a corpus, where the interesting shots are the ones whose
// `max_u Y(u)` straddles `T`: certified must be *exactly* "completed and within the horizon".
TEST(StockOnHCertificate, C1CertifiedIffCompletedWithinTheHorizon) {
    Corpus corpus = small_corpus(120);
    auto reference_mwpm = corpus.to_mwpm();
    size_t over_horizon_completions = 0;
    for (double t_edges : {0.5, 1.0, 1.5, 2.0}) {
        SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, t_edges);
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
        pm::total_weight_int horizon = decoder.horizon;

        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs);
            ASSERT_EQ(actual.weight, expected.weight);

            bool completed = profile.h_no_perfect_matching == 0;
            bool within = profile.max_dual_at_completion <= horizon;
            ASSERT_EQ(profile.certified == 1, completed && within)
                << "invariant 3: certified is not exactly (completed and max Y <= T) at T = " << t_edges;
            ASSERT_EQ(profile.escalated, profile.certified == 0);
            if (completed && !within)
                over_horizon_completions++;
        }
    }
    // If this is ever zero the sweep stopped exercising the case invariant 4 exists for, and the
    // test above is carrying it alone.
    std::printf(
        "[ INFO     ] %zu completing-but-over-T shots seen across the horizon sweep\n", over_horizon_completions);
}

// ---------------------------------------------------------------------------------------------
// C2 — no perfect matching in `H` is reported, not thrown, and the shot escalates to stock's exact
// answer. `H` having no perfect matching is a valid escalation trigger (§M7.0), not a failure.
TEST(StockOnHCertificate, C2NoPerfectMatchingIsReportedAndEscalates) {
    pm::UserGraph graph = fixture_graph();
    pm::Mwpm reference = graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false);
    auto decoder = BallDecoder::from_mwpm(graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false), fixture_config());

    // `c` and `d` are 10 apart, the filter admits 8, and the graph has no boundary at all — so `H`
    // is two isolated nodes.
    std::vector<uint64_t> shot = {2, 3};
    std::vector<CommittedPair> committed;
    BallProfile profile;
    Phase1Outcome outcome;
    ASSERT_NO_THROW(outcome = decoder.decode_phase1_production_to_match_edges(shot, committed, &profile));
    ASSERT_EQ(profile.h_no_perfect_matching, 1);
    ASSERT_EQ(profile.certified, 0);
    ASSERT_EQ(outcome.status, TimelineStatus::TRUNCATED);
    // The certificate is undefined rather than merely unmeasured when `H` cannot complete.
    ASSERT_EQ(profile.max_dual_at_completion, 0);

    Answer expected = stock_decode(reference, shot, reference.flooder.graph.num_observables);
    ASSERT_EQ(expected.weight, in_time_units(reference, 10.0));
}

// The decoder-level statement of the same thing: the shot escalates and the output is stock's.
TEST(StockOnHCertificate, C2NoPerfectMatchingEndToEnd) {
    pm::UserGraph graph = fixture_graph();
    pm::Mwpm reference = graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false);
    BallConfig ball_config = fixture_config();
    auto decoder = BallDecoder::from_mwpm(graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false), ball_config);

    for (const std::vector<uint64_t>& shot : {std::vector<uint64_t>{2, 3}, std::vector<uint64_t>{0, 1, 2, 3}}) {
        BallProfile profile;
        Phase1Outcome outcome = decoder.decode_phase1_production(shot, &profile);
        ASSERT_EQ(outcome.status, TimelineStatus::TRUNCATED);
        // Whatever the trigger, the escalating branch leaves nothing behind: the next shot on the
        // same instance has to be as clean as the first.
        ASSERT_TRUE(decoder.h_mwpm.mwpm.node_arena.live.empty());
        ASSERT_TRUE(decoder.h_mwpm.mwpm.flooder.region_arena.live.empty());
    }
    (void)reference;
}

// ---------------------------------------------------------------------------------------------
// C4 — `T = infinity` and over-large `T`. An unbounded horizon is unsupported on `H`: with a finite
// `R` the tables cannot supply the edges the run would need. It is routed through the
// `phase1_on_ball_graph = false` oracle path, and `T <= T_max` fires before anything else.
TEST(StockOnHCertificate, C4UnboundedAndOversizedHorizonsAreRejected) {
    Corpus corpus = small_corpus(1);
    auto reference_mwpm = corpus.to_mwpm();

    SpecMatchingConfig unbounded = stock_config_for(reference_mwpm.flooder.graph, 1.0);
    unbounded.unbounded_horizon = true;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, unbounded, NUM_DISTINCT_WEIGHTS), std::invalid_argument);

    // Stock blossom with no horizon is only defined on `G`, where it is simply the inherited
    // decoder; asking for it on `H` is a configuration error rather than a silent fallback.
    SpecMatchingConfig on_g = stock_config_for(reference_mwpm.flooder.graph, 1.0);
    on_g.phase1_on_ball_graph = false;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, on_g, NUM_DISTINCT_WEIGHTS), std::invalid_argument);

    SpecMatchingConfig oversized = stock_config_for(reference_mwpm.flooder.graph, 1.0);
    oversized.T = 2 * oversized.ball.T_max;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, oversized, NUM_DISTINCT_WEIGHTS), std::invalid_argument);

    // The M1 harvest oracle has nothing to reproduce on this path, so asking for it is refused
    // rather than silently compared against a harvest that never runs (§M7's verification note).
    SpecMatchingConfig verified = stock_config_for(reference_mwpm.flooder.graph, 1.0);
    verified.verify_against_g = true;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, verified, NUM_DISTINCT_WEIGHTS), std::invalid_argument);
}

// The two harvest-flavoured ball entry points are closed on this path, because their callers read
// the escalation decision off an empty residual — which is exactly wrong for a completing-but-over-T
// shot.
TEST(StockOnHCertificate, HarvestEntryPointsAreClosedOnTheStockPath) {
    pm::UserGraph graph = fixture_graph();
    auto decoder = BallDecoder::from_mwpm(graph.to_mwpm(NUM_DISTINCT_WEIGHTS, false), fixture_config());
    std::vector<uint64_t> shot = {0, 1, 2, 3};
    std::vector<CommittedPair> committed;
    ASSERT_THROW(decoder.decode_phase1(shot), std::invalid_argument);
    ASSERT_THROW(decoder.decode_phase1_to_match_edges(shot, committed), std::invalid_argument);
}

// ---------------------------------------------------------------------------------------------
// C5 — no production horizon. `GraphFlooder::horizon` holds its disabled sentinel on every shot, on
// both graphs, and — in a debug build, where the gate is instrumented — the horizon gate is never
// consulted at all. The second half is the stronger statement: it says no *code* gated, not merely
// that the field ended up back at its sentinel.
TEST(StockOnHCertificate, C5NoProductionHorizon) {
    Corpus corpus = small_corpus(60);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, 1.5);
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

#ifndef NDEBUG
    pm::horizon_gate_stats.clear();
#endif
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        spec_matching_decode(decoder, shot, &profile);
        ASSERT_EQ(decoder.ball->h_mwpm.mwpm.flooder.horizon, pm::NO_HORIZON)
            << "invariant 1: a finite horizon reached the stock-on-H flooder";
        ASSERT_EQ(decoder.g_mwpm().flooder.horizon, pm::NO_HORIZON);
    }
#ifndef NDEBUG
    ASSERT_EQ(pm::horizon_gate_stats.admitted, 0u) << "invariant 1: the horizon gate ran on the stock path";
    ASSERT_EQ(pm::horizon_gate_stats.rejected, 0u) << "invariant 1: the horizon gate rejected an event";
#endif
}

// ---------------------------------------------------------------------------------------------
// C6 — SCAN and BITSET build the same `H`, so they must reach the same certified/escalate decision
// and the same output on every shot.
TEST(StockOnHCertificate, C6ScanAndBitsetAgree) {
    Corpus corpus = small_corpus(60);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig scan_config = stock_config_for(reference_mwpm.flooder.graph, 1.5);
    SpecMatchingConfig bitset_config = scan_config;
    bitset_config.mode = BallGraphBuildMode::BITSET;
    auto scan = SpecMatchingDecoder::from_detector_error_model(corpus.dem, scan_config, NUM_DISTINCT_WEIGHTS);
    auto bitset = SpecMatchingDecoder::from_detector_error_model(corpus.dem, bitset_config, NUM_DISTINCT_WEIGHTS);

    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile scan_profile;
        SpecMatchingProfile bitset_profile;
        Answer scan_answer = spec_matching_decode(scan, shot, &scan_profile);
        Answer bitset_answer = spec_matching_decode(bitset, shot, &bitset_profile);
        ASSERT_EQ(scan.ball->arena.graph.h_to_det, bitset.ball->arena.graph.h_to_det);
        ASSERT_EQ(scan_profile.certified, bitset_profile.certified);
        ASSERT_EQ(scan_profile.h_no_perfect_matching, bitset_profile.h_no_perfect_matching);
        ASSERT_EQ(scan_profile.max_dual_at_completion, bitset_profile.max_dual_at_completion);
        ASSERT_EQ(scan_answer.obs, bitset_answer.obs);
        ASSERT_EQ(scan_answer.weight, bitset_answer.weight);
    }
}

// ---------------------------------------------------------------------------------------------
// C7 — determinism and reusability. Two runs are bit-identical, and a shot decoded after a certified
// shot and after an escalated shot both match a fresh decoder exactly. Both `Mwpm` instances are
// left clean.
TEST(StockOnHCertificate, C7DeterminismAndReuse) {
    Corpus corpus = small_corpus(80);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, 1.0);
    auto first = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
    auto second = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    // A long-lived decoder that has seen everything, against a fresh one that has seen nothing.
    std::vector<Answer> streamed;
    std::vector<bool> escalated;
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        streamed.push_back(spec_matching_decode(first, shot, &profile));
        escalated.push_back(profile.escalated);
        ASSERT_TRUE(first.ball->h_mwpm.mwpm.node_arena.live.empty());
        ASSERT_TRUE(first.ball->h_mwpm.mwpm.flooder.region_arena.live.empty());
        ASSERT_TRUE(first.g_mwpm().node_arena.live.empty());
    }
    for (size_t i = 0; i < corpus.shots.size(); i++) {
        SpecMatchingProfile profile;
        Answer again = spec_matching_decode(second, corpus.shots[i], &profile);
        ASSERT_EQ(again.obs, streamed[i].obs) << "shot " << i;
        ASSERT_EQ(again.weight, streamed[i].weight) << "shot " << i;
        ASSERT_EQ(profile.escalated, escalated[i]) << "shot " << i;
    }

    // The pointed version of the same claim: whichever branch the previous shot took, the next shot
    // is decided as if the decoder were fresh.
    size_t after_certified = 0;
    size_t after_escalated = 0;
    for (size_t i = 1; i < corpus.shots.size(); i++) {
        auto fresh = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
        Answer alone = spec_matching_decode(fresh, corpus.shots[i]);
        ASSERT_EQ(alone.obs, streamed[i].obs) << "shot " << i << " decoded differently in a stream";
        ASSERT_EQ(alone.weight, streamed[i].weight) << "shot " << i;
        (escalated[i - 1] ? after_escalated : after_certified)++;
        if (after_escalated > 0 && after_certified > 0 && i > 20)
            break;  // Both predecessors exercised; the rest is the loop above's job.
    }
    ASSERT_GT(after_certified, 0u);
}

// ---------------------------------------------------------------------------------------------
// C8 — negative-weight DEMs. The preamble is exercised, and escalation is handed the **raw**
// syndrome rather than `H`'s post-preamble effective set (invariant 14) — which is the thing a
// negative-weight fixture, and only a negative-weight fixture, distinguishes.
TEST(StockOnHCertificate, C8NegativeWeightDems) {
    Corpus corpus = load_surface_code_d13_negative_weights(60);
    auto reference_mwpm = corpus.to_mwpm();
    ASSERT_FALSE(reference_mwpm.flooder.negative_weight_detection_events.empty())
        << "this corpus is supposed to have negative weight edges";

    for (double t_edges : {0.75, 2.0}) {
        SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, t_edges);
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
        size_t escalated = 0;
        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << "T = " << t_edges;
            ASSERT_EQ(actual.weight, expected.weight) << "T = " << t_edges;
            escalated += profile.escalated ? 1 : 0;
        }
        std::printf(
            "[ INFO     ] negative-weight corpus at T = %.2f edges: %zu/%zu escalated\n",
            t_edges,
            escalated,
            corpus.shots.size());
    }
}

// ---------------------------------------------------------------------------------------------
// C9 — the `q = 0` path. At a horizon where nothing escalates, no `Mwpm(G)` work happens at all:
// `escalation_ns` is identically zero and every shot is certified.
TEST(StockOnHCertificate, C9ZeroEscalationPathCostsNothing) {
    Corpus corpus = small_corpus(60);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, 6.0, 6.0);
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        Answer actual = spec_matching_decode(decoder, shot, &profile);
        Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
        ASSERT_EQ(actual.obs, expected.obs);
        ASSERT_EQ(actual.weight, expected.weight);
        ASSERT_EQ(profile.certified, 1) << "a shot escalated at T = 6 edge weights on a d = 5 corpus";
        ASSERT_FALSE(profile.escalated);
        ASSERT_EQ(profile.escalation_ns, 0);
        ASSERT_EQ(profile.stock_ns, 0);
    }
    ASSERT_EQ(decoder.stats.shots_escalated, 0u);
}

// ---------------------------------------------------------------------------------------------
// C10 — fuzz. Random horizons below `T_max`, random syndromes, random small DEMs. Level 1 and the
// escalation-trigger invariant have to hold in all of them, not only on the corpora the design was
// written against.
TEST(StockOnHCertificate, C10Fuzz) {
    std::mt19937_64 rng(20260907);
    for (size_t trial = 0; trial < 6; trial++) {
        size_t distance = 3 + 2 * (trial % 3);
        double noise = 0.005 + 0.005 * (double)(trial % 3);
        Corpus corpus = generate_surface_code_corpus(distance, 3, noise, 25, 4242 + trial);
        auto reference_mwpm = corpus.to_mwpm();
        double unit = edge_weight_units(reference_mwpm.flooder.graph);

        double t_max_edges = 3.0;
        std::uniform_real_distribution<double> horizon_dist(0.2, t_max_edges);
        SpecMatchingConfig config;
        config.T = horizon_dist(rng) * unit;
        config.ball.T_max = t_max_edges * unit;
        config.ball.R = 2 * config.ball.T_max;
        config.stock_on_h = true;
        config.mode = (trial % 2) ? BallGraphBuildMode::BITSET : BallGraphBuildMode::SCAN;
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

        // Half the syndromes are the sampler's; the other half are random subsets of the detectors,
        // which reach parities and densities a physical sampler does not.
        std::vector<std::vector<uint64_t>> shots = corpus.shots;
        size_t num_detectors = reference_mwpm.flooder.graph.nodes.size();
        std::uniform_int_distribution<size_t> count_dist(0, std::min<size_t>(12, num_detectors));
        std::uniform_int_distribution<size_t> det_dist(0, num_detectors - 1);
        for (size_t i = 0; i < 25; i++) {
            std::vector<uint64_t> shot;
            size_t n = count_dist(rng);
            for (size_t k = 0; k < n; k++)
                shot.push_back(det_dist(rng));
            std::sort(shot.begin(), shot.end());
            shot.erase(std::unique(shot.begin(), shot.end()), shot.end());
            shots.push_back(shot);
        }

        for (const auto& shot : shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << "trial " << trial << ", T = " << config.T;
            ASSERT_EQ(actual.weight, expected.weight) << "trial " << trial << ", T = " << config.T;
            bool completed = profile.h_no_perfect_matching == 0;
            ASSERT_EQ(profile.certified == 1, completed && profile.max_dual_at_completion <= decoder.horizon);
            ASSERT_EQ(profile.escalated, profile.certified == 0);
        }
    }
}

// ---------------------------------------------------------------------------------------------
// §M7.0's corollary, measured rather than asserted: every shot this scheme escalates is one the
// landed truncated scheme also escalates, so `q_this <= q_current` on identical shots.
TEST(StockOnHCertificate, EscalatesNoMoreOftenThanTheTruncatedScheme) {
    Corpus corpus = small_corpus(120);
    auto reference_mwpm = corpus.to_mwpm();
    for (double t_edges : {1.0, 1.5, 2.0}) {
        SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, t_edges);
        config.measure_truncated_reference = true;
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            spec_matching_decode(decoder, shot, &profile);
            decoder.stats.accumulate(profile);
            ASSERT_TRUE(profile.truncated_reference_measured);
            // The corollary, per shot: `max_u Y(u) <= cur_time` throughout, so a shot that completes
            // within `T` also completed *by* `T`. Certified therefore implies the truncated scheme
            // kept it too.
            if (profile.certified)
                ASSERT_FALSE(profile.truncated_reference_escalates)
                    << "a certified shot that the truncated scheme escalated contradicts §M7.0's corollary";
        }
        SpecMatchingSummary summary = summarize(decoder.stats);
        ASSERT_EQ(summary.shots_with_truncated_reference, corpus.shots.size());
        ASSERT_LE(summary.q, summary.q_current_on_same_corpus);
        std::printf(
            "[ INFO     ] T = %.2f edges: q_this = %.4f, q_current = %.4f over %llu shots\n",
            t_edges,
            summary.q,
            summary.q_current_on_same_corpus,
            (unsigned long long)summary.shots);
    }
}

// ---------------------------------------------------------------------------------------------
// The dual scan itself: `max_nested_dual` and the `blossom_parent` chain walk it is checked against
// agree, and the value is bounded by the clock — which is the property §M7.0's corollary rests on.
TEST(StockOnHCertificate, MaxDualAgreesWithTheChainWalk) {
    Corpus corpus = small_corpus(40);
    pm::Mwpm mwpm = corpus.to_mwpm();
    for (const auto& shot : corpus.shots) {
        if (pm::process_timeline_until_completion_or_report(mwpm, shot) != pm::CompletionStatus::COMPLETE) {
            mwpm.reset();
            continue;
        }
        pm::total_weight_int scanned = max_nested_dual(mwpm, shot);
        pm::total_weight_int walked = max_nested_dual_by_chain_walk(mwpm, shot);
        ASSERT_EQ(scanned, walked);
        ASSERT_LE(scanned, (pm::total_weight_int)mwpm.flooder.queue.cur_time)
            << "max_u Y(u) outran the clock, which would void the terminal certificate";
        mwpm.reset();
    }
}

// ---------------------------------------------------------------------------------------------
// The edges flavour on this path: the lifted correction's outer syndrome is the input syndrome, and
// its weight is the matching weight, on certified and escalated shots alike.
TEST(StockOnHCertificate, EdgesFlavor) {
    Corpus corpus = small_corpus(50);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = stock_config_for(reference_mwpm.flooder.graph, 1.5);
    config.edges_flavor = true;
    config.ball.store_paths = true;
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    std::vector<int64_t> edges;
    pm::total_weight_int weight = 0;
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        decoder.decode_to_edges(shot, edges, weight, &profile);

        std::map<uint64_t, int> parity;
        for (size_t i = 0; i < edges.size() / 2; i++) {
            for (int64_t node : {edges[2 * i], edges[2 * i + 1]}) {
                if (node >= 0)
                    parity[(uint64_t)node] ^= 1;
            }
        }
        std::vector<uint64_t> produced;
        for (const auto& entry : parity) {
            if (entry.second)
                produced.push_back(entry.first);
        }
        std::vector<uint64_t> expected = shot;
        std::sort(expected.begin(), expected.end());
        ASSERT_EQ(produced, expected) << "the lifted correction does not produce the input syndrome";

        Answer obs_answer = stock_decode(reference_mwpm, shot, decoder.num_observables);
        ASSERT_EQ(weight, obs_answer.weight);
    }
}
