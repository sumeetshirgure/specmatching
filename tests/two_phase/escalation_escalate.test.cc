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

/// §M3.3 X1–X9: stock escalation, and §M3.4's extract-only production path.
///
/// The whole correctness argument for M3 is one assertion — on an escalating shot the decoder's
/// output is bit-exact against stock on the same syndrome — so most of this file is about arranging
/// for escalating shots to exist and then making that assertion in as many configurations as the
/// design enumerates.

#include <cstdio>
#include <random>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/driver/two_phase_decoding.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

TwoPhaseConfig config_for(const pm::MatchingGraph& graph, double t_edges, double t_max_edges = 2.0) {
    double unit = edge_weight_units(graph);
    TwoPhaseConfig config;
    config.T = t_edges * unit;
    config.ball.T_max = t_max_edges * unit;
    config.ball.R = 2.0 * t_max_edges * unit;
    return config;
}

struct Answer {
    std::vector<uint8_t> obs;
    pm::total_weight_int weight{0};

    bool operator==(const Answer& rhs) const {
        return obs == rhs.obs && weight == rhs.weight;
    }
};

Answer stock_decode(pm::Mwpm& mwpm, const std::vector<uint64_t>& dets, size_t num_observables) {
    Answer answer;
    answer.obs.assign(num_observables, 0);
    pm::decode_detection_events(mwpm, dets, answer.obs.data(), answer.weight, false);
    return answer;
}

Answer two_phase_decode(TwoPhaseDecoder& decoder, const std::vector<uint64_t>& dets, TwoPhaseProfile* prof = nullptr) {
    Answer answer;
    answer.obs.assign(decoder.num_observables, 0);
    decoder.decode_to_obs(dets, answer.obs.data(), answer.weight, prof);
    return answer;
}

/// A corpus small enough that compiling ball tables per test is affordable, and dense enough that
/// something escalates at a sane horizon.
Corpus small_corpus(size_t shots = 60) {
    return generate_surface_code_corpus(5, 5, 0.01, shots, 20260804);
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// X5 — forced escalation. Written first, per the implementation order: with `T` small enough that
// every shot truncates, the decoder is bit-identical to stock on the entire corpus. This is the
// cheap end-to-end proof of X1, and it exercises the escalation path on every single shot.
TEST(EscalationEscalate, X5ForcedEscalationIsIdenticalToStock) {
    Corpus corpus = small_corpus();
    auto reference_mwpm = corpus.to_mwpm();
    // A quarter of an edge weight: no region can reach a partner, so every shot with any defect
    // leaves a surviving tree.
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 0.25);
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t escalated = 0;
    size_t non_empty = 0;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        Answer actual = two_phase_decode(decoder, shot, &profile);
        Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
        ASSERT_EQ(actual.obs, expected.obs) << "escalated shot diverged from stock in the observables";
        ASSERT_EQ(actual.weight, expected.weight) << "escalated shot diverged from stock in the weight";
        escalated += profile.escalated ? 1 : 0;
        non_empty += shot.empty() ? 0 : 1;
    }
    ASSERT_GT(non_empty, 0u);
    ASSERT_EQ(escalated, non_empty) << "at T = 0.25 edge weights every shot with a defect must escalate";
}

// ---------------------------------------------------------------------------------------------
// X1 — escalation identity, the deliverable. At the operating horizon only a few shots escalate,
// so this asserts the same equality as X5 *and* checks that escalating shots actually occurred.
TEST(EscalationEscalate, X1EscalationIdentityAtTheOperatingHorizon) {
    Corpus corpus = generate_surface_code_corpus(7, 7, 0.02, 120, 20260805);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t escalated = 0;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        Answer actual = two_phase_decode(decoder, shot, &profile);
        Answer expected = stock_decode(reference_mwpm, shot, decoder.num_observables);
        // Exact on *every* shot, escalating or not: that is what M3 buys.
        ASSERT_EQ(actual.obs, expected.obs);
        ASSERT_EQ(actual.weight, expected.weight);
        escalated += profile.escalated ? 1 : 0;
    }
    std::printf(
        "[ INFO     ] X1: %zu of %zu shots escalated at T = 1.0 edge weights\n", escalated, corpus.shots.size());
    ASSERT_GT(escalated, 0u) << "this corpus was chosen to produce escalating shots; it produced none";
}

// ---------------------------------------------------------------------------------------------
// X2 — trigger exactness. Escalation fires iff `residual.size() > 0`, asserted against the harvest
// itself rather than against a recomputed predicate. Needs the verification path, which is the only
// configuration that builds a residual on a truncated shot at all.
TEST(EscalationEscalate, X2TriggerFiresExactlyOnANonEmptyResidual) {
    Corpus corpus = small_corpus();
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 0.75);
    config.full_harvest_for_verification = true;
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        two_phase_decode(decoder, shot, &profile);
        ASSERT_EQ(profile.escalated, profile.residual_size > 0) << "invariant 12";
        ASSERT_EQ(profile.escalated, profile.truncated) << "invariant 12";
        ASSERT_EQ(profile.residual_size, profile.num_trees) << "one exposed defect per surviving tree";
        if (!profile.escalated)
            ASSERT_EQ(profile.escalation_ns, 0);
    }
}

// ---------------------------------------------------------------------------------------------
// X3 — full-syndrome input. Stock is handed the shot's raw detection events, not the residual and
// not §M2.1's post-preamble effective set. The negative-weight DEM is the meaningful case, because
// it is the one where the two differ: if the driver passed the seeded set, the answer would be
// stock's answer to a *different* syndrome, which is a well-formed wrong answer.
TEST(EscalationEscalate, X3StockSeesTheRawSyndrome) {
    Corpus corpus = load_surface_code_d13_negative_weights(40);
    auto reference_mwpm = corpus.to_mwpm();
    auto seeded_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 0.25);
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
    ASSERT_FALSE(decoder.g_mwpm().flooder.negative_weight_detection_events.empty())
        << "this corpus was chosen for its negative weights";

    size_t shots_where_the_two_differ = 0;
    std::vector<uint64_t> seeded;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        Answer actual = two_phase_decode(decoder, shot, &profile);
        ASSERT_TRUE(profile.escalated);
        Answer on_raw = stock_decode(reference_mwpm, shot, decoder.num_observables);
        ASSERT_EQ(actual.obs, on_raw.obs);
        ASSERT_EQ(actual.weight, on_raw.weight);

        decoder.ball->compute_seeded_detection_events(shot, seeded);
        std::vector<uint64_t> sorted_shot = shot;
        std::sort(sorted_shot.begin(), sorted_shot.end());
        if (seeded == sorted_shot)
            continue;
        shots_where_the_two_differ++;
        // The two inputs really do give different answers, so agreeing with the raw one is
        // evidence and not a coincidence.
        Answer on_seeded = stock_decode(seeded_mwpm, seeded, decoder.num_observables);
        if (!(on_seeded == on_raw))
            ASSERT_FALSE(actual == on_seeded) << "stock was handed the post-preamble set, not the raw syndrome";
    }
    ASSERT_GT(shots_where_the_two_differ, 0u) << "the preamble never changed the event set; the test proved nothing";
}

// ---------------------------------------------------------------------------------------------
// X4 — instance reusability, in both directions. Both `Mwpm` instances have to come out of an
// escalation clean: the one on `H`, which was abandoned mid-solve, and the one on `G`, which ran a
// full stock decode.
TEST(EscalationEscalate, X4InstanceReusabilityAcrossEscalation) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.02, 80, 20260806);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    auto reused = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t transitions = 0;
    bool previous_escalated = false;
    for (size_t i = 0; i < corpus.shots.size(); i++) {
        TwoPhaseProfile profile;
        Answer actual = two_phase_decode(reused, corpus.shots[i], &profile);

        // A decoder that has seen nothing else, on the same shot.
        auto fresh = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
        Answer expected = two_phase_decode(fresh, corpus.shots[i]);
        ASSERT_EQ(actual.obs, expected.obs) << "shot " << i << " depended on the shots before it";
        ASSERT_EQ(actual.weight, expected.weight) << "shot " << i << " depended on the shots before it";

        if (i > 0 && profile.escalated != previous_escalated)
            transitions++;
        previous_escalated = profile.escalated;
    }
    ASSERT_GT(transitions, 0u) << "the corpus never crossed between escalating and non-escalating shots";
}

// ---------------------------------------------------------------------------------------------
// X6 — the `q = 0` path. At a horizon where nothing escalates the common case must not pay for the
// fallback existing: no escalation time, and no work on the `Mwpm` over `G` at all. The second half
// is asserted structurally — `G`'s region arena has never handed out a region — rather than by
// timing, which would not distinguish "cheap" from "absent".
TEST(EscalationEscalate, X6TheCommonCaseDoesNotTouchG) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.002, 60, 20260807);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 2.0);
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t escalated = 0;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        two_phase_decode(decoder, shot, &profile);
        escalated += profile.escalated ? 1 : 0;
        if (!profile.escalated)
            ASSERT_EQ(profile.escalation_ns, 0);
    }
    ASSERT_EQ(escalated, 0u) << "this operating point was chosen so that nothing escalates";
    ASSERT_TRUE(decoder.g_mwpm().flooder.region_arena.allocated.empty())
        << "the empty-residual path allocated a region on G; it must not touch G at all";
    ASSERT_TRUE(decoder.g_mwpm().node_arena.allocated.empty());
    ASSERT_EQ(decoder.g_mwpm().flooder.queue.cur_time, 0);
}

// ---------------------------------------------------------------------------------------------
// X7 — determinism, across runs and across a serialise/load of the ball artifact. An escalating
// shot goes through two solvers and a teardown, which is the most state a shot can touch.
TEST(EscalationEscalate, X7DeterminismAcrossRunsAndArtifactRoundTrip) {
    Corpus corpus = small_corpus(40);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 1.0);

    auto first = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
    std::string path = std::tmpnam(nullptr);
    first.save_ball_artifact(path);
    auto second = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);
    auto from_artifact =
        TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS, path.c_str());
    std::remove(path.c_str());

    size_t escalated = 0;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile profile;
        Answer a = two_phase_decode(first, shot, &profile);
        Answer b = two_phase_decode(second, shot);
        Answer c = two_phase_decode(from_artifact, shot);
        ASSERT_EQ(a.obs, b.obs);
        ASSERT_EQ(a.weight, b.weight);
        ASSERT_EQ(a.obs, c.obs) << "loading the ball artifact changed the answer";
        ASSERT_EQ(a.weight, c.weight);
        escalated += profile.escalated ? 1 : 0;
    }
    ASSERT_GT(escalated, 0u);
}

// ---------------------------------------------------------------------------------------------
// X8 — bypass identity. The extract-only production path and the full §M1.3 harvest must agree,
// bit for bit, on every non-escalating shot of the corpus, in both flavours. This is what licenses
// §M3.4: without it the bypass is an untested claim that two jobs compute nothing anyone reads.
TEST(EscalationEscalate, X8ExtractOnlyMatchesFullHarvest) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.01, 80, 20260808);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig production = config_for(reference_mwpm.flooder.graph, 1.0);
    production.edges_flavor = true;
    production.ball.store_paths = true;
    TwoPhaseConfig verification = production;
    verification.full_harvest_for_verification = true;

    auto fast = TwoPhaseDecoder::from_detector_error_model(corpus.dem, production, NUM_DISTINCT_WEIGHTS);
    auto slow = TwoPhaseDecoder::from_detector_error_model(corpus.dem, verification, NUM_DISTINCT_WEIGHTS);

    size_t compared = 0;
    for (const auto& shot : corpus.shots) {
        TwoPhaseProfile fast_profile;
        TwoPhaseProfile slow_profile;
        Answer fast_answer = two_phase_decode(fast, shot, &fast_profile);
        Answer slow_answer = two_phase_decode(slow, shot, &slow_profile);
        ASSERT_EQ(fast_profile.escalated, slow_profile.escalated) << "the two paths disagreed on the branch itself";
        if (fast_profile.escalated)
            continue;
        compared++;
        ASSERT_EQ(fast_answer.obs, slow_answer.obs) << "obs flavour: extract-only diverged from full harvest";
        ASSERT_EQ(fast_answer.weight, slow_answer.weight) << "obs flavour: weights differ";
        ASSERT_EQ(fast_profile.dual_sum_at_truncation, slow_profile.dual_sum_at_truncation)
            << "extract-only dropped the dual sum the §M4.2 certificate reads";
        ASSERT_EQ(fast_profile.committed_pairs_frozen, slow_profile.committed_pairs_frozen);
        ASSERT_EQ(fast_profile.committed_boundary, slow_profile.committed_boundary);

        std::vector<int64_t> fast_edges;
        std::vector<int64_t> slow_edges;
        pm::total_weight_int fast_weight = 0;
        pm::total_weight_int slow_weight = 0;
        fast.decode_to_edges(shot, fast_edges, fast_weight);
        slow.decode_to_edges(shot, slow_edges, slow_weight);
        ASSERT_EQ(fast_edges, slow_edges) << "edges flavour: extract-only diverged from full harvest";
        ASSERT_EQ(fast_weight, slow_weight);
    }
    ASSERT_GT(compared, 0u);
}

// ---------------------------------------------------------------------------------------------
// X9 — the branch is O(1) and correct.
//
// Correctness is asserted inside `run_timeline` on every shot of every test in this binary: the
// `live` vector's `empty()`, the free-list comparison, and a full sweep over live regions all have
// to agree. What is asserted here is the other half — that the production path really does perform
// no per-tree work to reach the branch, counted rather than timed, so no amount of scheduler noise
// can make it pass by accident.
TEST(EscalationEscalate, X9ProductionPathDoesNoTreeWork) {
    Corpus corpus = generate_surface_code_corpus(5, 5, 0.01, 80, 20260809);
    auto reference_mwpm = corpus.to_mwpm();
    TwoPhaseConfig config = config_for(reference_mwpm.flooder.graph, 1.0);
    auto decoder = TwoPhaseDecoder::from_detector_error_model(corpus.dem, config, NUM_DISTINCT_WEIGHTS);

    size_t escalated = 0;
    size_t completed = 0;
    for (const auto& shot : corpus.shots) {
        decoder.ball->harvester.counters.reset();
        TwoPhaseProfile profile;
        two_phase_decode(decoder, shot, &profile);
        const HarvestCounters& counters = decoder.ball->harvester.counters;

        ASSERT_EQ(counters.full_harvests, 0u) << "the production path ran the full harvest";
        ASSERT_EQ(counters.tree_nodes_visited, 0u) << "the production path enumerated alternating tree nodes";
        ASSERT_EQ(counters.base_descents, 0u) << "M1.4's base descent is still on the production path";
        if (profile.escalated) {
            escalated++;
            ASSERT_EQ(counters.extract_only_harvests, 0u) << "an escalating shot harvested instead of being abandoned";
            ASSERT_EQ(counters.extractions, 0u);
        } else {
            completed++;
            // One extract-only harvest, unless §A's resolver settled the whole of `H` off the
            // solver — then §A.4 skips the build, the solve and the extraction together, and zero
            // harvests is the honest count rather than a missed one. Stated as an equality against
            // which case the shot was in, so a harvest that goes missing for any *other* reason
            // still fails here.
            size_t h_nodes = decoder.ball->arena.graph.num_nodes();
            bool fully_resolved = (size_t)profile.defects_resolved_small == h_nodes;
            ASSERT_EQ(counters.extract_only_harvests, fully_resolved ? 0u : 1u)
                << "defects_resolved_small = " << profile.defects_resolved_small << " of " << h_nodes << " H nodes";
        }
    }
    ASSERT_GT(escalated, 0u);
    ASSERT_GT(completed, 0u);

    // And the same question asked the slow way, on a hand-driven shot, so that the equality is
    // checked in a release build too rather than only under `assert`.
    auto& h_mwpm = decoder.ball->h_mwpm.mwpm;
    ASSERT_EQ(any_alternating_tree_survives(h_mwpm), any_alternating_tree_survives_by_sweep(h_mwpm));
}
