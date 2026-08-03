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

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_TRUNCATION_TRUNCATED_TIMELINE_H
