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

#include "specmatching/spec_matching/certificate/max_dual.h"

#include <algorithm>
#include <cassert>

#include "specmatching/sparse_blossom/driver/mwpm_decoding.h"
#include "specmatching/sparse_blossom/flooder/graph.h"
#include "specmatching/spec_matching/perf/spec_matching_profile.h"
#include "specmatching/spec_matching/truncation/exposed_blossom.h"
#include "specmatching/spec_matching/truncation/harvest.h"

namespace pm {
namespace spec_matching {

namespace {

/// Folds `f` over every detection event the shot's dual solution covers.
///
/// The graph's own negative-weight detection events are included alongside the caller's, exactly as
/// the stock extraction path treats them: they are real detection events with real regions. On `H`
/// that list is empty by construction (§M2.1 runs the preamble on `G` first and `H` never sees a
/// negative weight), so this costs nothing there and keeps the routine correct if it is ever
/// pointed at `G`.
template <typename Fn>
void for_each_dual_carrying_node(const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, const Fn& f) {
    const auto& nodes = mwpm.flooder.graph.nodes;
    const std::vector<uint64_t>* event_lists[2] = {&detection_events, &mwpm.flooder.negative_weight_detection_events};
    for (const std::vector<uint64_t>* events : event_lists) {
        for (uint64_t det : *events) {
            if (det >= nodes.size())
                continue;
            const pm::DetectorNode& node = nodes[det];
            if (node.region_that_arrived == nullptr)
                continue;
            f(node);
        }
    }
}

}  // namespace

pm::total_weight_int max_nested_dual(const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    const pm::cumulative_time_int time = mwpm.flooder.queue.cur_time;
    pm::total_weight_int max_dual = 0;
    for_each_dual_carrying_node(mwpm, detection_events, [&](const pm::DetectorNode& node) {
        max_dual = std::max(max_dual, (pm::total_weight_int)nested_dual_sum(node, time));
    });
    return max_dual;
}

pm::total_weight_int max_nested_dual_by_chain_walk(
    const pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events) {
    const pm::cumulative_time_int time = mwpm.flooder.queue.cur_time;
    pm::total_weight_int max_dual = 0;
    for_each_dual_carrying_node(mwpm, detection_events, [&](const pm::DetectorNode& node) {
        // `Y(u) = Sum_{S ∋ u} y_S`, one term per level of blossom nesting. A blossom's `radius` is
        // the *extra* growth since it was formed, not a total, so the chain sums rather than
        // maximises — and every region below the top is frozen, so evaluating each term at `time`
        // and adding is the same number as adding the `Varying`s and evaluating once.
        pm::total_weight_int dual = 0;
        for (const pm::GraphFillRegion* region = node.region_that_arrived; region != nullptr;
             region = region->blossom_parent) {
            dual += (pm::total_weight_int)region->radius.get_distance_at_time(time);
        }
        max_dual = std::max(max_dual, dual);
    });
    return max_dual;
}

CertificateOutcome run_stock_and_certify(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, horizon_int horizon, bool time_dual_scan) {
    // Debug invariant 1, on the way in: the production flooder must never observe a finite horizon.
    // This is the whole point of M7 — the flooder is stock, so a finite horizon here is a bug.
    assert(mwpm.flooder.horizon == pm::NO_HORIZON && "invariant 1: a horizon leaked into the stock-on-H flooder");

    CertificateOutcome outcome;
    if (pm::process_timeline_until_completion_or_report(mwpm, detection_events) ==
        pm::CompletionStatus::NO_PERFECT_MATCHING) {
        // Correct, and expected: no perfect matching in `H` implies some dual of the true optimum
        // exceeds `T` (§M7.0), so the shot escalates. The instance is left standing for the caller
        // to tear down.
        outcome.status = CertificateStatus::NO_PERFECT_MATCHING;
        assert(mwpm.flooder.horizon == pm::NO_HORIZON);
        return outcome;
    }

    // Debug invariant 5: a completed run has no surviving alternating tree, asked both the O(1) way
    // and the slow way, so the two cannot drift apart unnoticed.
    assert(!any_alternating_tree_survives(mwpm) && !any_alternating_tree_survives_by_sweep(mwpm));

    // The one terminal scan. Once, after the loop reaches `NO_EVENT`, over the completed state and
    // before anything is extracted — the nesting this reads is what extraction destroys.
    HiResTimer timer;
    if (time_dual_scan)
        timer.start();
    outcome.max_dual = max_nested_dual(mwpm, detection_events);
    if (time_dual_scan)
        outcome.dual_scan_ns = timer.elapsed_ns();

    // Debug invariant 2: the same maximum by an independent route.
    assert(
        outcome.max_dual == max_nested_dual_by_chain_walk(mwpm, detection_events) &&
        "invariant 2: the dual readout and the blossom_parent chain walk disagree on max_u Y(u)");

    outcome.status = outcome.max_dual <= (pm::total_weight_int)horizon ? CertificateStatus::CERTIFIED
                                                                       : CertificateStatus::DUAL_EXCEEDS_HORIZON;
    assert(mwpm.flooder.horizon == pm::NO_HORIZON);
    return outcome;
}

}  // namespace spec_matching
}  // namespace pm
