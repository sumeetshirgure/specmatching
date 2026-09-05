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

/// §2 of the component profiler — the one clock `component_profiler` reads, and the only timing
/// code in `benchmarks/spec_matching/component_profiler.cc`.
///
/// Two entry points, both free functions rather than a stopwatch object, because §2 asks for
/// **two reads per component** and a `now` is the shape that makes that countable at the call site:
///
///   - `hires_now_ns()`   — this thread's run time so far, in nanoseconds;
///   - `hires_timer_name()` — which clock produced it, written verbatim to `run.log`.
///
/// Plus `hires_timer_overhead_ns()`, §2's startup calibration.
///
/// ## Platform
///
/// Detected at compile time, and an unsupported platform is a build error rather than a silent
/// fallback to a wall clock — a wall-clock fallback does not fail, it just quietly charges every
/// component for whatever else the machine was doing, which is the one failure mode here that
/// produces plausible numbers instead of obviously broken ones.
///
///   - **macOS**: `clock_gettime(CLOCK_THREAD_CPUTIME_ID)`, this thread's user + kernel
///     nanoseconds, served by the `thread_selfusage` mach trap. The value derives from the 24 MHz
///     `cntvct_el0` counter, so granularity is ~42 ns and a read costs a few hundred ns.
///     `hires_timer_name()` is `thread_cputime_macos`.
///   - **Linux**: the latency profiler's existing Linux path, unchanged — `perf/thread_timer.h`'s
///     backend selection, which prefers the per-thread `perf_event` cycle counter read with
///     `rdpmc` (~8 ns/read) and degrades to `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` where the
///     kernel refuses the event. `hires_timer_name()` is whatever that header settled on, so a
///     `run.log` records which of the two produced its rows.
///
/// Both platforms land on a **thread-scoped** clock: it advances only while this thread is on a
/// CPU, so a component the scheduler interrupts is not charged the scheduler's time. That is what
/// §2's "nothing further is done about preemption" rests on. `hires_timer_is_thread_scoped()` says
/// whether the claim actually holds on this run — on Linux the backend is chosen at run time and a
/// machine that refuses both thread-scoped clocks would land on wall time, which `run.log` has to
/// record rather than assume away.

#ifndef SPECMATCHING_BENCHMARKS_SPEC_MATCHING_HIRES_TIMER_H
#define SPECMATCHING_BENCHMARKS_SPEC_MATCHING_HIRES_TIMER_H

#include <cstdint>

#if defined(__APPLE__)
#include <time.h>
#elif defined(__linux__)
#include "specmatching/perf/thread_timer.h"
#else
#error "unsupported platform: component_profiler needs a thread-scoped nanosecond clock"
#endif

namespace pm {
namespace spec_matching {
namespace profiler {

/// This thread's run time so far, in nanoseconds. Only differences are meaningful; the origin is
/// unspecified and differs by backend.
inline uint64_t hires_now_ns() {
#if defined(__APPLE__)
    struct timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
#else
    // The Linux path, taken through `thread_timer.h` rather than reimplemented: whichever backend
    // that header selected, its ticks scaled by its own calibration. `ThreadTimer` itself is a
    // start/elapsed pair and this needs a bare `now` — for the calibration below as much as for the
    // per-component reads — so it reads the same two pieces `ThreadTimer` does.
    //
    // On the cycle-counter backend a tick is a cycle and `ns_per_tick` is a fraction, so the cast
    // truncates: a delta is off by at most 1 ns, against a read that costs ~8. On the `cputime`
    // backend ticks are already nanoseconds and `ns_per_tick` is exactly 1.
    const pm::perf::internal::ThreadClock& clock = pm::perf::internal::thread_clock();
    return (uint64_t)((double)pm::perf::internal::read_ticks(clock) * clock.ns_per_tick);
#endif
}

/// The clock's name, for `run.log`.
inline const char* hires_timer_name() {
#if defined(__APPLE__)
    return "thread_cputime_macos";
#else
    return pm::perf::current_backend_name();
#endif
}

/// Does the clock stop while this thread is off the CPU? Always true on macOS, where there is one
/// backend and it is thread-scoped; on Linux it is a property of what the run-time backend
/// selection found, so it is reported rather than assumed.
inline bool hires_timer_is_thread_scoped() {
#if defined(__APPLE__)
    return true;
#else
    return pm::perf::backend_is_thread_scoped(pm::perf::current_backend());
#endif
}

/// §2's overhead calibration: 10^6 back-to-back `hires_now_ns()` calls, `(last - first) / 10^6`.
///
/// Run once at startup and written to `run.log` as `timer_overhead_ns`. Rows are stored
/// **uncorrected** — the analysis subtracts this — because a correction baked into the rows cannot
/// be undone by a reader who disagrees with it.
///
/// Returned as a double: on the Linux cycle-counter backend the true figure is single-digit
/// nanoseconds and an integer division would quantise it to uselessness.
inline double hires_timer_overhead_ns() {
    constexpr uint64_t CALLS = 1000000;
    uint64_t first = hires_now_ns();
    // Volatile so the loop is 10^6 real readings rather than something the optimiser folds into
    // one; the clock's own cost is exactly what is being measured.
    volatile uint64_t last = first;
    for (uint64_t i = 0; i < CALLS; i++)
        last = hires_now_ns();
    return (double)(last - first) / (double)CALLS;
}

}  // namespace profiler
}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_BENCHMARKS_SPEC_MATCHING_HIRES_TIMER_H
