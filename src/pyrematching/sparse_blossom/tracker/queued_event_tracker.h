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

#ifndef PYREMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H
#define PYREMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H

#include <limits>

#include "pyrematching/sparse_blossom/tracker/radix_heap_queue.h"

namespace pm {

/// Sentinel horizon meaning "no truncation": every event is admitted into the queue.
/// A `GraphFlooder` carries this value except while `process_timeline_until_horizon` is running,
/// so the stock decode path never observes a finite horizon.
constexpr cumulative_time_int NO_HORIZON = std::numeric_limits<cumulative_time_int>::max();

/// Debug-only instrumentation for the M1 horizon gate.
///
/// Only updated in debug builds (`NDEBUG` undefined). It exists to discharge the debug-build
/// invariant "nothing with cumulative time > horizon is ever enqueued": after a truncated shot,
/// `rejected` counts the gated-out schedulings and `max_admitted_time` must be <= the horizon.
/// Only *gated* call sites (those that pass a horizon) are counted; shrink events are exempt from
/// the gate by design and are therefore not counted here.
struct HorizonGateStats {
    uint64_t admitted{0};
    uint64_t rejected{0};
    cumulative_time_int max_admitted_time{std::numeric_limits<cumulative_time_int>::min()};

    void clear() {
        admitted = 0;
        rejected = 0;
        max_admitted_time = std::numeric_limits<cumulative_time_int>::min();
    }
};

/// Process-wide debug counters for the horizon gate. Not thread safe; debug builds only.
inline HorizonGateStats horizon_gate_stats{};

/// This class is responsible for ensuring that a "look at me!" event is in the event queue.
///
/// This object also attempts to avoid spamming the event queue with redundant events. If there are
/// two reasons to look at an object at a specific time, it will only put one event into the queue.
/// If an object is already going to be looked at at time T, and it wants to be looked at at time
/// T+2, the look-event for time T+2 will not be enqueued right away. Instead, as part of processing
/// the time T event, the tracker will take care of enqueueing the time T+2 event.
struct QueuedEventTracker {
    cyclic_time_int desired_time{0};
    cyclic_time_int queued_time{0};
    bool has_desired_time{false};
    bool has_queued_time{false};

    /// Resets the tracker to its initial idle state.
    void clear() {
        desired_time = cyclic_time_int{0};
        queued_time = cyclic_time_int{0};
        has_desired_time = {false};
        has_queued_time = {false};
    }

    /// Tells the tracker a desired look-at-me event. The tracker will handle inserting an event
    /// into the event queue, if necessary.
    ///
    /// `horizon` truncates the timeline (M1): an event that would occur strictly after the horizon
    /// is *never inserted into the queue* — not inserted and later skipped — and the tracker is
    /// left with no desired event, exactly as if the object had no upcoming event at all. The
    /// default is `NO_HORIZON`, under which this is bit-for-bit the stock behaviour.
    ///
    /// `cumulative_event_time` is the un-wrapped time of `ev`, and is what the horizon is compared
    /// against. `ev.time` is a `cyclic_time_int` and cannot be ordered against a cumulative bound,
    /// so callers must pass the cumulative time they computed *before* wrapping it.
    template <bool use_validation>
    inline void set_desired_event(
        FloodCheckEvent ev,
        radix_heap_queue<use_validation> &queue,
        cumulative_time_int cumulative_event_time = 0,
        cumulative_time_int horizon = NO_HORIZON) {
        if (cumulative_event_time > horizon) {
            has_desired_time = false;
#ifndef NDEBUG
            horizon_gate_stats.rejected++;
#endif
            return;
        }
#ifndef NDEBUG
        if (horizon != NO_HORIZON) {
            horizon_gate_stats.admitted++;
            horizon_gate_stats.max_admitted_time =
                std::max(horizon_gate_stats.max_admitted_time, cumulative_event_time);
        }
#endif
        has_desired_time = true;
        desired_time = ev.time;
        if (!has_queued_time || queued_time > ev.time) {
            queued_time = ev.time;
            has_queued_time = true;
            queue.enqueue(ev);
        }
    }

    /// Indicates that it's no longer necessary to look at the object at a later time.
    inline void set_no_desired_event() {
        has_desired_time = false;
    }

    /// Notifies the tracker that a relevant look-at-me event has been dequeued from the event
    /// queue. The result of this method is whether or not to discard the event (due to it being
    /// no longer desired) instead of continuing processing it. Returning false means discard,
    /// returning true means keep. This method also handles requeueing another look-at-me event
    /// if doing so was deferred while the earlier event was in the queue.
    template <bool use_validation>
    inline bool dequeue_decision(FloodCheckEvent ev, radix_heap_queue<use_validation> &queue) {
        // Only the most recent event this tracker put into the queue is valid. Older events
        // are forgotten and must not be processed, because otherwise an event storm can be
        // created as stale events trigger redundant enqueues.
        if (!has_queued_time || ev.time != queued_time) {
            return false;
        }
        has_queued_time = false;

        // If the event isn't for the CURRENTLY desired look-at-me time, discard it.
        if (!has_desired_time) {
            return false;
        }
        if (ev.time != desired_time) {
            // Requeue the event if the desired time is a little later.
            has_queued_time = true;
            queued_time = desired_time;
            ev.time = desired_time;
            queue.enqueue(ev);
            return false;
        }

        // All systems go! Process away!
        has_desired_time = false;
        return true;
    }
};

}  // namespace pm

#endif  // PYREMATCHING_FILL_MATCH_QUEUED_EVENT_TRACKER_H
