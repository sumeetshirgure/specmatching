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

/// §M4.2 and M5: the end-to-end driver, in both flavours.
///
/// The headline test is `ExactnessEndToEnd`, and it is an *equality*, not a ratio. The old plan
/// gated M4 on "LER vs exact, no catastrophic regression"; with stock escalation the output is
/// exact MWPM on every shot, so the gate is bit-equality against stock and a failure is a Phase-1
/// commit-policy bug, not an accuracy tradeoff.

#include <cstdio>
#include <map>
#include <random>

#include "gtest/gtest.h"

#include "specmatching/spec_matching/driver/spec_matching_decoding.h"
#include "tests/spec_matching/spec_matching_test_util.h"

using namespace pm::spec_matching;
using namespace pm::spec_matching::test;

namespace {

SpecMatchingConfig config_for(const pm::MatchingGraph& graph, double t_edges, double t_max_edges = 2.0) {
    double unit = edge_weight_units(graph);
    SpecMatchingConfig config;
    config.T = t_edges * unit;
    config.ball.T_max = t_max_edges * unit;
    config.ball.R = 2.0 * t_max_edges * unit;
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

std::vector<int64_t> sorted_edges(std::vector<int64_t> edges) {
    std::vector<std::pair<int64_t, int64_t>> pairs;
    for (size_t i = 0; i < edges.size() / 2; i++) {
        int64_t u = edges[2 * i];
        int64_t v = edges[2 * i + 1];
        if (v >= 0 && v < u)
            std::swap(u, v);
        pairs.emplace_back(u, v);
    }
    std::sort(pairs.begin(), pairs.end());
    edges.clear();
    for (const auto& pair : pairs) {
        edges.push_back(pair.first);
        edges.push_back(pair.second);
    }
    return edges;
}

/// The syndrome an edge set produces: every detector an odd number of the edges touch. A boundary
/// half-edge touches one detector. This is `H @ edges` computed independently of the decoder, which
/// is what makes it evidence that the correction is valid rather than merely well-formed.
std::vector<uint64_t> syndrome_of_edges(const std::vector<int64_t>& edges) {
    std::map<uint64_t, int> parity;
    for (size_t i = 0; i < edges.size() / 2; i++) {
        for (int64_t node : {edges[2 * i], edges[2 * i + 1]}) {
            if (node >= 0)
                parity[(uint64_t)node] ^= 1;
        }
    }
    std::vector<uint64_t> out;
    for (const auto& entry : parity) {
        if (entry.second)
            out.push_back(entry.first);
    }
    return out;
}

/// The observables an edge set crosses, read off the search graph rather than off anything the
/// decoder produced.
std::vector<uint8_t> observables_of_edges(
    const pm::SearchGraph& graph, const std::vector<int64_t>& edges, size_t num_observables) {
    std::vector<uint8_t> obs(num_observables, 0);
    for (size_t i = 0; i < edges.size() / 2; i++) {
        const pm::SearchDetectorNode& u = graph.nodes[(size_t)edges[2 * i]];
        int64_t v = edges[2 * i + 1];
        size_t k = v < 0 ? u.index_of_neighbor(nullptr)
                         : u.index_of_neighbor(const_cast<pm::SearchDetectorNode*>(&graph.nodes[(size_t)v]));
        for (size_t o : u.neighbor_observable_indices[k])
            obs[o] ^= 1;
    }
    return obs;
}

Corpus edges_corpus(size_t shots = 60) {
    return generate_surface_code_corpus(5, 5, 0.01, shots, 20260810);
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// `T = infinity` identity, through the driver, on the oracle front end. This is M1 test 1 restated
// one layer up: the whole pipeline at an unbounded horizon has to be the stock decoder.
TEST(SpecMatchingDriver, UnboundedHorizonIsIdenticalToStock) {
    for (Corpus corpus : {load_surface_code_d13(60), load_surface_code_d13_negative_weights(60)}) {
        auto reference_mwpm = corpus.to_mwpm();
        SpecMatchingConfig config;
        config.phase1_on_ball_graph = false;
        config.unbounded_horizon = true;
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << corpus.name;
            ASSERT_EQ(actual.weight, expected.weight) << corpus.name;
            ASSERT_FALSE(profile.truncated) << "an unbounded horizon cannot truncate";
            ASSERT_FALSE(profile.escalated);
        }
    }
}

// `T = infinity` on the ball graph is rejected outright rather than left to fail invariant 6.
TEST(SpecMatchingDriver, UnboundedHorizonOnTheBallGraphIsRejected) {
    Corpus corpus = edges_corpus(1);
    SpecMatchingConfig config;
    config.phase1_on_ball_graph = true;
    config.unbounded_horizon = true;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS), std::invalid_argument);
}

// ---------------------------------------------------------------------------------------------
// Exactness end to end. Bit-exact against stock on **every** shot of the corpus, escalating or not,
// at the operating horizon. This replaces the old LER-ratio gate: the output is exact MWPM, so the
// gate is equality.
TEST(SpecMatchingDriver, ExactnessEndToEnd) {
    for (Corpus corpus :
         {load_surface_code_d13(60), load_surface_code_d13_negative_weights(60), load_toric_code_d5(60)}) {
        auto reference_mwpm = corpus.to_mwpm();
        SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 2.0);
        auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

        size_t escalated = 0;
        for (const auto& shot : corpus.shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << corpus.name
                                                << ": if this fails the bug is in Phase 1's commit "
                                                   "policy, not in the escalation";
            ASSERT_EQ(actual.weight, expected.weight) << corpus.name;
            escalated += profile.escalated ? 1 : 0;
        }
        std::printf(
            "[ INFO     ] %s: %zu of %zu shots escalated at T = 2 edge weights\n",
            corpus.name.c_str(),
            escalated,
            corpus.shots.size());
    }
}

// ---------------------------------------------------------------------------------------------
// The per-shot certificate. `dual_sum_at_truncation` is a lower bound on the exact optimum, so the
// ratio `weight_out / dual_sum` is at least 1 — a cheap independent check on the dual that costs
// one division. It is expected to be *tight* on non-escalating shots, which is what the reported
// distribution is for.
TEST(SpecMatchingDriver, PerShotCertificate) {
    Corpus corpus = load_surface_code_d13(80);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 2.0);
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    double worst = 0;
    double sum = 0;
    size_t counted = 0;
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        spec_matching_decode(decoder, shot, &profile);
        if (profile.escalated) {
            // The truncated dual went with the rest of Phase 1's discarded result; there is nothing
            // to certify against, and reporting the discarded number would invite a meaningless
            // ratio.
            ASSERT_EQ(profile.dual_sum_at_truncation, 0);
            continue;
        }
        if (profile.dual_sum_at_truncation == 0)
            continue;
        double ratio = (double)profile.weight_out / (double)profile.dual_sum_at_truncation;
        ASSERT_GE(profile.weight_out, profile.dual_sum_at_truncation)
            << "the committed weight fell below Sum_S y_S, which is a lower bound on the optimum";
        worst = std::max(worst, ratio);
        sum += ratio;
        counted++;
    }
    ASSERT_GT(counted, 0u);
    std::printf(
        "[ INFO     ] certificate ratio over %zu non-escalating shots: mean %.4f, max %.4f\n",
        counted,
        sum / (double)counted,
        worst);
}

// ---------------------------------------------------------------------------------------------
// Batch determinism, and that `decode_batch` agrees with the per-shot entry point it wraps.
TEST(SpecMatchingDriver, BatchDeterminism) {
    Corpus corpus = edges_corpus(60);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    auto batched = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
    auto per_shot = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t num_obs = batched.num_observables;
    std::vector<uint8_t> first(corpus.shots.size() * num_obs, 0);
    std::vector<pm::total_weight_int> first_weights(corpus.shots.size(), 0);
    batched.decode_batch(corpus.shots, first.data(), first_weights.data(), true);

    std::vector<uint8_t> second(corpus.shots.size() * num_obs, 0);
    std::vector<pm::total_weight_int> second_weights(corpus.shots.size(), 0);
    batched.decode_batch(corpus.shots, second.data(), second_weights.data(), true);
    ASSERT_EQ(first, second);
    ASSERT_EQ(first_weights, second_weights);

    for (size_t i = 0; i < corpus.shots.size(); i++) {
        Answer expected = spec_matching_decode(per_shot, corpus.shots[i]);
        std::vector<uint8_t> actual(first.begin() + (long)(i * num_obs), first.begin() + (long)((i + 1) * num_obs));
        ASSERT_EQ(actual, expected.obs) << "shot " << i;
        ASSERT_EQ(first_weights[i], expected.weight) << "shot " << i;
    }

    // The aggregate accumulator saw both passes and the summary reconciles against it.
    ASSERT_EQ(batched.stats.shots, 2 * corpus.shots.size());
    SpecMatchingSummary summary = summarize(batched.stats);
    ASSERT_EQ(summary.shots, batched.stats.shots);
    ASSERT_DOUBLE_EQ(summary.q, (double)batched.stats.shots_escalated / (double)batched.stats.shots);
}

// ---------------------------------------------------------------------------------------------
// M5. The lifted edge set is a valid correction: its outer syndrome equals the input syndrome, over
// the full corpus rather than a sample.
TEST(SpecMatchingDriver, EdgesFlavorOuterSyndromeMatchesTheInput) {
    Corpus corpus = edges_corpus(80);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    config.edges_flavor = true;
    config.ball.store_paths = true;
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t escalated = 0;
    std::vector<int64_t> edges;
    pm::total_weight_int weight = 0;
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        decoder.decode_to_edges(shot, edges, weight, &profile);
        std::vector<uint64_t> expected = shot;
        std::sort(expected.begin(), expected.end());
        ASSERT_EQ(syndrome_of_edges(edges), expected) << "the emitted edge set is not a correction for this syndrome";
        escalated += profile.escalated ? 1 : 0;
    }
    ASSERT_GT(escalated, 0u) << "the escalating branch of the edges flavour was never exercised";
}

// Observables derived from the edge set match the obs-flavour output on the same shot. The two go
// through completely different readouts — ball-table observable id lists against ball-table stored
// paths — so agreeing is a real check and not a tautology.
TEST(SpecMatchingDriver, EdgesFlavorObservablesMatchTheObsFlavor) {
    Corpus corpus = edges_corpus(80);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    config.edges_flavor = true;
    config.ball.store_paths = true;
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    std::vector<int64_t> edges;
    pm::total_weight_int edge_weight = 0;
    for (const auto& shot : corpus.shots) {
        Answer obs_answer = spec_matching_decode(decoder, shot);
        decoder.decode_to_edges(shot, edges, edge_weight, nullptr);
        std::vector<uint8_t> from_edges =
            observables_of_edges(decoder.g_mwpm().search_flooder.graph, edges, decoder.num_observables);
        ASSERT_EQ(from_edges, obs_answer.obs) << "the two flavours disagree on the correction's homology";
        ASSERT_EQ(edge_weight, obs_answer.weight) << "the two flavours disagree on the matching weight";
    }
}

// Edge-flavour identity at `T = infinity`, through the oracle front end.
TEST(SpecMatchingDriver, EdgesFlavorUnboundedHorizonIsIdenticalToStock) {
    Corpus corpus = edges_corpus(60);
    auto reference_mwpm = corpus.to_mwpm(true);
    SpecMatchingConfig config;
    config.phase1_on_ball_graph = false;
    config.unbounded_horizon = true;
    config.edges_flavor = true;
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    auto obs_reference_mwpm = corpus.to_mwpm();
    std::vector<int64_t> actual;
    std::vector<int64_t> expected;
    pm::total_weight_int weight = 0;
    for (const auto& shot : corpus.shots) {
        decoder.decode_to_edges(shot, actual, weight, nullptr);
        expected.clear();
        pm::decode_detection_events_to_edges(reference_mwpm, shot, expected);
        ASSERT_EQ(sorted_edges(actual), sorted_edges(expected));
        ASSERT_EQ(weight, stock_decode(obs_reference_mwpm, shot, decoder.num_observables).weight);
    }
}

// Edge-flavour escalation identity: on every escalating shot the edge list is bit-exact against
// stock (§M3.3 X1, edges flavour). Forced, so it is every shot.
TEST(SpecMatchingDriver, EdgesFlavorEscalationIsIdenticalToStock) {
    Corpus corpus = edges_corpus(60);
    auto reference_mwpm = corpus.to_mwpm(true);
    SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, 0.25);
    config.edges_flavor = true;
    config.ball.store_paths = true;
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    auto obs_reference_mwpm = corpus.to_mwpm();
    std::vector<int64_t> actual;
    std::vector<int64_t> expected;
    pm::total_weight_int weight = 0;
    size_t escalated = 0;
    for (const auto& shot : corpus.shots) {
        SpecMatchingProfile profile;
        decoder.decode_to_edges(shot, actual, weight, &profile);
        if (!profile.escalated)
            continue;
        escalated++;
        expected.clear();
        pm::decode_detection_events_to_edges(reference_mwpm, shot, expected);
        ASSERT_EQ(sorted_edges(actual), sorted_edges(expected));
        // The weight an escalated shot reports is the *matching* weight, which is stock's obs-flavour
        // weight. Reading it off the emitted edge set instead would differ whenever two matched
        // pairs' paths share an edge and cancel it out of the correction.
        Answer obs_answer = stock_decode(obs_reference_mwpm, shot, decoder.num_observables);
        ASSERT_EQ(weight, obs_answer.weight);
    }
    ASSERT_GT(escalated, 0u);
}

// ---------------------------------------------------------------------------------------------
// More than 64 observables. §0 requires every new path to work for an arbitrary observable count,
// and above 64 the `obs_int` mask is unusable — so the committed pairs are read out through the
// ball tables' observable *id* lists instead. That is a different route to the same answer and it
// has no other test, because every corpus in `data/` is a single-observable memory experiment.
TEST(SpecMatchingDriver, WideObservableCount) {
    std::string text = "error(0.05) D0 L0\n";
    const size_t chain = 70;
    for (size_t i = 0; i + 1 < chain; i++) {
        text += "error(0.1) D" + std::to_string(i) + " D" + std::to_string(i + 1) + " L" + std::to_string(i + 1) + "\n";
    }
    text += "error(0.05) D" + std::to_string(chain - 1) + " L" + std::to_string(chain) + "\n";
    stim::DetectorErrorModel dem(text.c_str());
    ASSERT_GT(dem.count_observables(), sizeof(pm::obs_int) * 8);

    auto reference_mwpm = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, true);

    // Syndromes are sampled by flipping random *errors* and taking their endpoints, not by picking
    // random detectors: on a chain, independent detectors are far apart and every shot would
    // escalate, leaving the committed readout — the thing this test exists for — untried.
    std::mt19937_64 rng(20260811);
    std::vector<std::vector<uint64_t>> shots;
    for (size_t trial = 0; trial < 200; trial++) {
        std::vector<int> parity(chain, 0);
        for (size_t k = 0; k < 3; k++) {
            size_t error_index = rng() % (chain + 1);
            if (error_index > 0)
                parity[error_index - 1] ^= 1;
            if (error_index < chain)
                parity[error_index] ^= 1;
        }
        shots.emplace_back();
        for (size_t d = 0; d < chain; d++) {
            if (parity[d])
                shots.back().push_back(d);
        }
    }

    // Both readouts, deliberately: at `T = 2` edge weights every shot completes and the answer comes
    // out of the ball tables' observable id lists; at a quarter of an edge weight every shot
    // escalates and it comes out of stock. Both have to equal stock, and both are on the wide path.
    for (double t_edges : {2.0, 0.25}) {
        SpecMatchingConfig config = config_for(reference_mwpm.flooder.graph, t_edges);
        auto decoder = SpecMatchingDecoder::from_detector_error_model(dem, config, NUM_DISTINCT_WEIGHTS);
        size_t escalated = 0;
        for (const auto& shot : shots) {
            SpecMatchingProfile profile;
            Answer actual = spec_matching_decode(decoder, shot, &profile);
            Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
            ASSERT_EQ(actual.obs, expected.obs) << "wide-observable readout diverged from stock at T = " << t_edges;
            ASSERT_EQ(actual.weight, expected.weight) << "at T = " << t_edges;
            escalated += profile.escalated ? 1 : 0;
        }
        std::printf(
            "[ INFO     ] wide observables at T = %.2f: %zu of %zu shots escalated\n",
            t_edges,
            escalated,
            shots.size());
        if (t_edges > 1.0)
            ASSERT_LT(escalated, shots.size()) << "the wide-observable committed readout never ran";
        else
            ASSERT_GT(escalated, 0u) << "the wide-observable escalation path never ran";
    }
}

// The edges flavour refuses to run without the machinery it needs, rather than emitting an empty
// or half-lifted edge set.
TEST(SpecMatchingDriver, EdgesFlavorRequiresItsConfiguration) {
    Corpus corpus = edges_corpus(1);
    auto reference_mwpm = corpus.to_mwpm();
    SpecMatchingConfig no_paths = config_for(reference_mwpm.flooder.graph, 1.0);
    no_paths.edges_flavor = true;
    ASSERT_THROW(
        SpecMatchingDecoder::from_detector_error_model(corpus.dem, no_paths, NUM_DISTINCT_WEIGHTS), std::invalid_argument);

    SpecMatchingConfig obs_only = config_for(reference_mwpm.flooder.graph, 1.0);
    auto decoder = SpecMatchingDecoder::from_detector_error_model(corpus.dem, obs_only, NUM_DISTINCT_WEIGHTS);
    std::vector<int64_t> edges;
    pm::total_weight_int weight = 0;
    ASSERT_THROW(decoder.decode_to_edges(corpus.shots[0], edges, weight, nullptr), std::invalid_argument);
}
