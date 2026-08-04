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

/// §M2.9.7 H6 — the E1 and E2 experiments of §M2.9.5, as standalone reported measurements.
///
/// These run **before** any deferred-shatter code exists, which is the order the design requires:
/// the saving §M2.9.5 describes is gated on both experiments, and if E2 fails the shatter is
/// required even in the obs flavour.
///
/// The first test below is the one that makes the other two mean anything: it checks that the
/// read-only evaluator, driven with upstream's own split rule, reproduces what the real destructive
/// extraction produces, region for region. Without that, a "no difference" result could just mean
/// the evaluator is not modelling extraction.

#include <algorithm>
#include <cstdio>

#include "gtest/gtest.h"

#include "pyrematching/two_phase/truncation/harvest.h"
#include "pyrematching/two_phase/truncation/split_experiment.h"
#include "pyrematching/two_phase/truncation/truncated_timeline.h"
#include "tests/two_phase/two_phase_test_util.h"

using namespace pm::two_phase;
using namespace pm::two_phase::test;

namespace {

std::vector<Corpus> split_corpora() {
    std::vector<Corpus> corpora;
    corpora.push_back(load_surface_code_d13(40));
    corpora.push_back(load_surface_code_d13_negative_weights(40));
    corpora.push_back(load_toric_code_d5(40));
    corpora.push_back(generate_surface_code_corpus(9, 9, 0.01, 40, 7));
    return corpora;
}

/// The top-level matched regions of a truncated state, one entry per pair.
std::vector<pm::GraphFillRegion*> matched_representatives(pm::Mwpm& mwpm) {
    std::vector<pm::GraphFillRegion*> out;
    for (pm::GraphFillRegion* region : mwpm.flooder.region_arena.live) {
        if (region->blossom_parent != nullptr || region->alt_tree_node != nullptr)
            continue;
        if (region->match.region == nullptr)
            continue;
        if (region->match.edge.loc_from < region->match.edge.loc_to)
            out.push_back(region);
    }
    return out;
}

}  // namespace

// The evaluator is faithful: driven with `SplitRule::UPSTREAM` it reproduces
// `Mwpm::shatter_blossom_and_extract_matches` exactly, on every matched pair of every corpus shot.
TEST(TwoPhaseSplitExperiment, EvaluatorReproducesTheRealExtraction) {
    size_t pairs_checked = 0;
    size_t blossom_pairs_checked = 0;
    for (const Corpus& corpus : split_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3, 6}) {
                process_timeline_until_horizon(mwpm, shot, unit * multiple);

                // Predict first — the evaluator is read-only, so the state it saw is the state the
                // extraction below starts from.
                std::vector<pm::GraphFillRegion*> representatives = matched_representatives(mwpm);
                std::vector<SplitEvaluation> predicted;
                predicted.reserve(representatives.size());
                for (pm::GraphFillRegion* region : representatives) {
                    if (!region->blossom_children.empty() || !region->match.region->blossom_children.empty())
                        blossom_pairs_checked++;
                    predicted.push_back(evaluate_shatter(*region, SplitRule::UPSTREAM));
                }

                for (size_t i = 0; i < representatives.size(); i++) {
                    pm::MatchingResult actual = mwpm.shatter_blossom_and_extract_matches(representatives[i]);
                    ASSERT_EQ(actual.obs_mask, predicted[i].obs_mask)
                        << corpus.name << ": the evaluator predicted a different observable mask";
                    ASSERT_EQ(actual.weight, predicted[i].weight)
                        << corpus.name << ": the evaluator predicted a different weight";
                    pairs_checked++;
                }

                // The state is now half-extracted, so it cannot be harvested; drop it and start the
                // next shot from a fresh instance.
                mwpm = corpus.to_mwpm();
            }
        }
    }
    ASSERT_GT(pairs_checked, 0u);
    ASSERT_GT(blossom_pairs_checked, 0u) << "no matched blossom ever occurred, so the split rule was never exercised";
}

// H6 / E1 and E2. Reported, not asserted-into-a-pass: the design asks for these two numbers, and
// what they are decides whether §M2.9.5 is open or closed. The only hard assertion is that the
// experiment was not vacuous.
TEST(TwoPhaseSplitExperiment, H6ReportE1AndE2OverTheCorpus) {
    SplitExperimentStats overall;
    for (const Corpus& corpus : split_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        SplitExperimentStats per_corpus;
        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3, 6}) {
                process_timeline_until_horizon(mwpm, shot, unit * multiple);
                probe_split_dependence(mwpm, per_corpus);
                harvester.harvest_to_obs(mwpm, shot);
            }
        }

        std::printf(
            "[ INFO     ] %s: %llu shots, %llu with a matched blossom, %llu splittable matches, max nesting depth "
            "%d; E1 obs moved %llu, interior obs non-zero %llu, E2 weight moved %llu, cycles with non-zero obs "
            "%llu/%llu\n",
            corpus.name.c_str(),
            (unsigned long long)per_corpus.shots,
            (unsigned long long)per_corpus.shots_with_matched_blossoms,
            (unsigned long long)per_corpus.splittable_matches,
            per_corpus.max_nesting_depth,
            (unsigned long long)per_corpus.obs_differed,
            (unsigned long long)per_corpus.interior_obs_nonzero,
            (unsigned long long)per_corpus.weight_differed,
            (unsigned long long)per_corpus.cycles_with_nonzero_obs,
            (unsigned long long)per_corpus.cycles_examined);
        overall.merge(per_corpus);
    }

    std::printf(
        "[ INFO     ] E1/E2 overall: %llu splittable matches over %llu shots; E1 failures %llu (%.4f), interior obs "
        "non-zero %llu (%.4f), E2 failures %llu (%.4f), cycles with non-zero obs %llu/%llu, max nesting depth %d\n",
        (unsigned long long)overall.splittable_matches,
        (unsigned long long)overall.shots,
        (unsigned long long)overall.obs_differed,
        overall.splittable_matches ? (double)overall.obs_differed / (double)overall.splittable_matches : 0.0,
        (unsigned long long)overall.interior_obs_nonzero,
        overall.splittable_matches ? (double)overall.interior_obs_nonzero / (double)overall.splittable_matches : 0.0,
        (unsigned long long)overall.weight_differed,
        overall.splittable_matches ? (double)overall.weight_differed / (double)overall.splittable_matches : 0.0,
        (unsigned long long)overall.cycles_with_nonzero_obs,
        (unsigned long long)overall.cycles_examined,
        overall.max_nesting_depth);

    ASSERT_GT(overall.splittable_matches, 0u) << "no matched blossom ever occurred, so E1/E2 measured nothing";
}

// E2, stated structurally rather than statistically, because the structural statement is the one
// that generalises past the corpus: the weight `shatter_blossom_and_extract_matches` accumulates is
// the sum of `y_S` over every region in the nesting, each counted exactly once, and no choice of
// split point changes which regions exist. Any counterexample here is a bug in the evaluator.
TEST(TwoPhaseSplitExperiment, E2CommittedWeightIsIndependentOfTheSplitPoint) {
    size_t compared = 0;
    for (const Corpus& corpus : split_corpora()) {
        auto mwpm = corpus.to_mwpm();
        Harvester harvester;
        auto unit = median_edge_weight(mwpm.flooder.graph);

        for (const auto& shot : corpus.shots) {
            for (int multiple : {1, 2, 3, 6}) {
                process_timeline_until_horizon(mwpm, shot, unit * multiple);
                for (pm::GraphFillRegion* region : matched_representatives(mwpm)) {
                    if (region->blossom_children.empty() && region->match.region->blossom_children.empty())
                        continue;
                    SplitEvaluation upstream = evaluate_shatter(*region, SplitRule::UPSTREAM);
                    SplitEvaluation rotated = evaluate_shatter(*region, SplitRule::ROTATED);
                    ASSERT_EQ(upstream.weight, rotated.weight)
                        << corpus.name << ": committed weight moved when the split point moved";
                    compared++;
                }
                harvester.harvest_to_obs(mwpm, shot);
            }
        }
    }
    ASSERT_GT(compared, 0u);
}
