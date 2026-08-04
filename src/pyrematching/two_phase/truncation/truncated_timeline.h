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

#ifndef PYREMATCHING_TWO_PHASE_TRUNCATION_TRUNCATED_TIMELINE_H
#define PYREMATCHING_TWO_PHASE_TRUNCATION_TRUNCATED_TIMELINE_H

#include <vector>

#include "pyrematching/sparse_blossom/matcher/mwpm.h"
#include "pyrematching/two_phase/truncation/horizon.h"

namespace pm {
namespace two_phase {

enum class TimelineStatus {
    /// The timeline ran to completion before the horizon: no alternating trees survive and the
    /// primal solution is a perfect matching. Harvesting yields an empty residual.
    COMPLETE,
    /// The timeline was truncated at the horizon with alternating trees still standing. This is an
    /// expected outcome, not an error.
    TRUNCATED,
};

/// Restores `flooder.horizon` to `pm::NO_HORIZON` on scope exit, including on an exception, so
/// that a throwing shot cannot leave a finite horizon behind for the next (stock) decode.
struct HorizonGuard {
    pm::GraphFlooder& flooder;

    HorizonGuard(pm::GraphFlooder& flooder, horizon_int horizon) : flooder(flooder) {
        flooder.horizon = horizon;
    }
    ~HorizonGuard() {
        flooder.horizon = pm::NO_HORIZON;
    }

    HorizonGuard(const HorizonGuard&) = delete;
    HorizonGuard& operator=(const HorizonGuard&) = delete;
};

/// Runs sparse blossom on `detection_events` until either the matching completes or the timeline
/// reaches `horizon`, whichever happens first.
///
/// Events beyond the horizon are never inserted into the event queue (see
/// `QueuedEventTracker::set_desired_event`), so truncation costs nothing beyond the gate itself;
/// the queue drains on its own.
///
/// Unlike `process_timeline_until_completion` this does *not* throw when alternating trees survive
/// — that is the `TRUNCATED` outcome, and the surviving trees are the input to `harvest`. It does
/// still propagate the argument errors of the shared preamble (empty-queue precondition,
/// out-of-range detection event index).
///
/// On return with `TRUNCATED`, `mwpm.flooder.queue.cur_time == horizon`: the duals have been
/// advanced to exactly `T`, which is what makes `Y(u) == T` hold for every residual defect.
/// The caller is expected to `harvest(...)` and then `mwpm.reset()`.
TimelineStatus process_timeline_until_horizon(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, horizon_int horizon);

/// §M2.9.6 measurement 4: the serial depth of the solve, in dependent events.
///
/// Harvest's critical path only means something next to the solve's, and "harvest is roughly a
/// third of what remains" was an estimate that the M2.9 exit checkpoint requires be replaced by a
/// measurement. Wall time will not do — it is dominated by this machine's cache hierarchy, and the
/// architecture the critical-path metric is *for* has a different one.
///
/// The model is a longest-chain count over the detector graph. Each detector node carries a depth;
/// an event that names nodes `u` and `v` has depth `1 + max(depth[u], depth[v])`, and writes that
/// back to both. A blossom shatter names no edge, so it takes the max over every node the blossom
/// owns and writes back to all of them. `depth` is then the length of the longest chain of events
/// that had to happen in order, which is what a machine with unbounded width would still pay.
struct TimelineDepthModel {
    /// Per-node chain depth, indexed by detector id. Reset per shot.
    std::vector<int> node_depth;
    /// Base of the graph's node vector, so that a `DetectorNode*` can be turned into an index.
    const pm::DetectorNode* node_base{nullptr};
    /// The shot's answer: the deepest chain of dependent events.
    int depth{0};
    /// Events processed, so that "depth vs count" — the available parallelism — can be reported.
    int events{0};

    void begin(const pm::MatchingGraph& graph);
    void observe(const pm::MwpmEvent& event);
};

/// As `process_timeline_until_horizon`, and additionally fills `depth_model`. Split out rather than
/// defaulted so that the measured path and the hot path are visibly the same loop with one extra
/// call in it.
TimelineStatus process_timeline_until_horizon_measured(
    pm::Mwpm& mwpm,
    const std::vector<uint64_t>& detection_events,
    horizon_int horizon,
    TimelineDepthModel& depth_model);

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_TRUNCATION_TRUNCATED_TIMELINE_H
