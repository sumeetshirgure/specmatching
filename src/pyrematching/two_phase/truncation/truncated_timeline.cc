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

#include "pyrematching/two_phase/truncation/truncated_timeline.h"

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"

namespace pm {
namespace two_phase {

TimelineStatus process_timeline_until_horizon(
    pm::Mwpm& mwpm, const std::vector<uint64_t>& detection_events, horizon_int horizon) {
    // Reuses the stock preamble verbatim: empty-queue precondition, cur_time = 0, and the
    // negative-weight detection event marking dance.
    pm::begin_timeline(mwpm, detection_events);

    // Exception-safe: the flooder goes back to the sentinel horizon however we leave this scope,
    // so the stock decode path never observes a finite horizon.
    HorizonGuard guard(mwpm.flooder, horizon);

    while (true) {
        auto event = mwpm.flooder.run_until_next_mwpm_notification();
        if (event.event_type == pm::NO_EVENT)
            break;
        mwpm.process_event(event);
    }

    // Alternating tree nodes are allocated from `node_arena` and returned to it when a tree is
    // resolved, so "some node is still checked out" is exactly "some alternating tree survives".
    bool trees_survive = mwpm.node_arena.allocated.size() != mwpm.node_arena.available.size();
    if (!trees_survive)
        return TimelineStatus::COMPLETE;

    if (horizon != pm::NO_HORIZON) {
        // Advance the clock to exactly `T`. The queue holds nothing at or before the horizon
        // (node events past it were never enqueued, and the loop above stopped at the first
        // exempt shrink event past it), so the surviving regions grow undisturbed to `T`. Doing
        // this here, once, is what makes `Y(u) == T` an exact equality at harvest.
        mwpm.flooder.queue.cur_time = horizon;
    }
    return TimelineStatus::TRUNCATED;
}

}  // namespace two_phase
}  // namespace pm
