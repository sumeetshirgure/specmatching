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

#ifndef PYREMATCHING_TWO_PHASE_TRUNCATION_SPLIT_EXPERIMENT_H
#define PYREMATCHING_TWO_PHASE_TRUNCATION_SPLIT_EXPERIMENT_H

#include <algorithm>
#include <cstdint>

#include "pyrematching/sparse_blossom/matcher/mwpm.h"

namespace pm {
namespace two_phase {

/// §M2.9.5, experiments E1 and E2. **Measurement only — there is no deferred-shatter code here, and
/// per the design there must not be until these two experiments have been reported.**
///
/// A `match` in sparse blossom is between *regions*, not defects. When a matched region is a blossom
/// the interior pairing is not stored anywhere: it is derived at extraction time by splitting the
/// tight odd cycle at the point where the external match attaches. The question §M2.9.5 is gated on
/// is whether the two things harvest actually reports — the committed observable flip and the
/// committed weight — depend on that split point. If neither does, `decode_to_obs` need not shatter
/// matched blossoms at all.
///
/// The machinery below answers it by *evaluating* an extraction rather than performing one: it
/// reproduces `Mwpm::shatter_blossom_and_extract_matches` exactly, as a pure function of the region
/// structure, and so can be run twice on the same live state with two different split rules. Nothing
/// it touches is mutated or freed.

/// Which member of each odd cycle is removed before the rest is paired up.
enum class SplitRule {
    /// Upstream's: the sub-blossom containing the node the external match attaches at.
    UPSTREAM,
    /// One step around the cycle from upstream's, at every level of nesting. This is the
    /// perturbation E1 and E2 are measured against.
    ROTATED,
};

/// What an extraction would contribute. The same two fields `pm::MatchingResult` carries, kept
/// separate so that it is obvious nothing here has been extracted for real.
struct SplitEvaluation {
    pm::obs_int obs_mask{0};
    pm::total_weight_int weight{0};

    bool operator==(const SplitEvaluation& rhs) const {
        return obs_mask == rhs.obs_mask && weight == rhs.weight;
    }
};

/// The (observable mask, weight) that shattering `region` — matched to `partner` via `edge`, or to
/// the boundary when `partner` is null — would produce under `rule`.
///
/// Mirrors `Mwpm::shatter_blossom_and_extract_matches` and
/// `Mwpm::pair_and_shatter_subblossoms_and_extract_matches` line for line, with the one substitution
/// that the split index comes from `rule` instead of always from the attach point. Read-only: the
/// region nesting, the matches and the arena are all left exactly as they were.
SplitEvaluation evaluate_shatter(
    const pm::GraphFillRegion& region,
    const pm::GraphFillRegion* partner,
    const pm::CompressedEdge& edge,
    SplitRule rule);

/// Convenience wrapper for a region that is already matched: reads the partner and edge off it.
SplitEvaluation evaluate_shatter(const pm::GraphFillRegion& region, SplitRule rule);

/// Tallies for E1 and E2, accumulated over a campaign.
struct SplitExperimentStats {
    uint64_t shots{0};
    /// Shots in which at least one matched blossom existed at truncation, i.e. shots where the
    /// question is not vacuous.
    uint64_t shots_with_matched_blossoms{0};
    /// Matched pairs at truncation in which at least one side is a blossom, so a split exists.
    uint64_t splittable_matches{0};
    /// E1, as the design words it: matches whose committed observable mask moved when the split
    /// point was rotated, with the external match edge held where it is.
    uint64_t obs_differed{0};
    /// E1 in the form that actually decides §M2.9.5: matches where the *interior* pairing
    /// contributes a non-zero observable flip, so that not shattering at all would give a different
    /// answer from shattering. This is the saving the section is after, stated directly.
    uint64_t interior_obs_nonzero{0};
    /// E2: matches whose committed weight moved when the split point was rotated.
    uint64_t weight_differed{0};
    /// Blossoms — at every level of nesting, on both sides of the match — whose tight odd cycle
    /// has a non-zero total observable mask. §M2.9.5's premise is that two splits differ by the
    /// full cycle and that such a cycle is homologically trivial far below the code distance, so
    /// this is the premise measured on its own.
    uint64_t cycles_examined{0};
    uint64_t cycles_with_nonzero_obs{0};
    /// Deepest nesting seen in a matched blossom, so that "the rotation was exercised at depth" is
    /// stated rather than assumed.
    int max_nesting_depth{0};

    void reset() {
        *this = SplitExperimentStats();
    }
    void merge(const SplitExperimentStats& other) {
        shots += other.shots;
        shots_with_matched_blossoms += other.shots_with_matched_blossoms;
        splittable_matches += other.splittable_matches;
        obs_differed += other.obs_differed;
        interior_obs_nonzero += other.interior_obs_nonzero;
        weight_differed += other.weight_differed;
        cycles_examined += other.cycles_examined;
        cycles_with_nonzero_obs += other.cycles_with_nonzero_obs;
        max_nesting_depth = std::max(max_nesting_depth, other.max_nesting_depth);
    }
};

/// Runs E1 and E2 over the state a truncated timeline left behind, and accumulates into `stats`.
///
/// Call it after `process_timeline_until_horizon` and **before** harvest: it reads the live regions
/// and does not mutate them, so harvesting afterwards produces exactly what it would have anyway.
/// Boundary matches and unmatched (tree) regions are skipped — a boundary match has no cycle to
/// split, and exposed root blossoms are explicitly outside E1/E2's scope (§M2.9.5's hard
/// exceptions).
void probe_split_dependence(pm::Mwpm& mwpm, SplitExperimentStats& stats);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_TRUNCATION_SPLIT_EXPERIMENT_H
