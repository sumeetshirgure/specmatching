// Copyright 2022 PyMatching Contributors
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

#ifndef PYREMATCHING2_MWPM_DECODING_H
#define PYREMATCHING2_MWPM_DECODING_H

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "stim.h"

namespace pm {

struct ExtendedMatchingResult {
    std::vector<uint8_t> obs_crossed;
    total_weight_int weight;
    ExtendedMatchingResult();
    explicit ExtendedMatchingResult(size_t num_observables);

    bool operator==(const ExtendedMatchingResult& rhs) const;

    bool operator!=(const ExtendedMatchingResult& rhs) const;

    ExtendedMatchingResult(std::vector<uint8_t> obs_crossed, total_weight_int weight);

    void reset();

    ExtendedMatchingResult& operator+=(const ExtendedMatchingResult& rhs);
    ExtendedMatchingResult operator+(const ExtendedMatchingResult& rhs) const;
};

inline void ExtendedMatchingResult::reset() {
    std::fill(obs_crossed.begin(), obs_crossed.end(), 0);
    weight = 0;
}

void fill_bit_vector_from_obs_mask(pm::obs_int obs_mask, uint8_t* obs_begin_ptr, size_t num_observables);
obs_int bit_vector_to_obs_mask(const std::vector<uint8_t>& bit_vector);

Mwpm detector_error_model_to_mwpm(
    const stim::DetectorErrorModel& detector_error_model,
    pm::weight_int num_distinct_weights,
    bool ensure_search_flooder_included = false,
    bool enable_correlations = false);

/// Seeds a fresh timeline: checks the queue is empty, resets `cur_time` to 0, and creates the
/// detection events for `detection_events`, including the marking/unmarking dance that cancels
/// detection events implied by negative weight edges.
///
/// Factored out of `process_timeline_until_completion` so that the truncated timeline
/// (`pm::spec_matching::process_timeline_until_horizon`) reuses this preamble verbatim rather than
/// copying it.
void begin_timeline(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Did the timeline complete, or is the graph one on which no perfect matching exists?
enum class CompletionStatus {
    /// The queue drained with no alternating tree standing: the primal is a perfect
    /// matching-with-boundary and the dual is optimal.
    COMPLETE,
    /// Alternating trees survived the drained queue, so the graph admits no perfect matching on
    /// this syndrome. An *expected* outcome for a caller that has somewhere to escalate to, which
    /// is why this is reported rather than thrown (§M7.2).
    NO_PERFECT_MATCHING,
};

/// Stock sparse blossom, run to completion, reporting rather than throwing.
///
/// Identical to `process_timeline_until_completion` in every respect except the failure mode: where
/// that one calls `Mwpm::reset` and throws `std::invalid_argument`, this returns
/// `NO_PERFECT_MATCHING` and **leaves the instance untouched**, so the caller owns the teardown and
/// can inspect the state first. `process_timeline_until_completion` is implemented on top of this,
/// so there is one timeline loop and one tree-survival test, not two.
///
/// No horizon is set, consulted or restored: `GraphFlooder::horizon` keeps its `pm::NO_HORIZON`
/// sentinel throughout, which is the machine-checkable form of §M7.1's "no hardcoded horizon".
CompletionStatus process_timeline_until_completion_or_report(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

MatchingResult decode_detection_events_for_up_to_64_observables(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, bool edge_correlations);

/// Used to decode detection events for an existing Mwpm object `mwpm', and a vector of
/// detection event indices `detection_events'. The predicted observables are XOR-ed into an
/// existing uint8_t array with at least `mwpm.flooder.graph.num_observables' elements,
/// the pointer to the first element of which is passed as the `obs_begin_ptr' argument.
/// The weight of the MWPM solution is added to the `weight' argument.
void decode_detection_events(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    uint8_t* obs_begin_ptr,
    pm::total_weight_int& weight,
    bool edge_correlations);

/// Decode detection events using a Mwpm object and vector of detection event indices
/// Returns the compressed edges in the matching: the pairs of detection events that are
/// matched to each other via paths.
void decode_detection_events_to_match_edges(pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events);

/// Turns `mwpm.flooder.match_edges` — matched detection-event *pairs* — into the *edges* of a
/// correction, appended to `edges` as detector id pairs with `-1` for the boundary. Flips each
/// edge along the search graph's shortest path between the matched pair, adds the graph's own
/// negative-weight edges, and drops any edge flipped an even number of times.
///
/// Factored out of `decode_detection_events_to_edges` (pyrematching M5) so that the spec-matching
/// driver's oracle front end can reuse it verbatim on a *truncated* timeline, where the match
/// edges come from `pm::spec_matching::harvest_to_match_edges` rather than from a completed decode.
/// Copying it instead would have put the cancellation pass in two places, and the cancellation is
/// the part with no independent check on it.
void expand_match_edges_to_edges(pm::Mwpm& mwpm, std::vector<int64_t>& edges);

/// Decode detection events using a Mwpm object and vector of detection event indices.
/// Returns the edges in the matching: these are pairs of *detectors* forming *edges* in the
/// matching solution (rather than pairs of detection *events* matched via *paths* as returned
/// instead by `decode_detection_events_to_match_edges`).
void decode_detection_events_to_edges(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);

void decode_detection_events_to_edges_with_edge_correlations(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, std::vector<int64_t>& edges);

}  // namespace pm

#endif  // PYREMATCHING2_MWPM_DECODING_H
