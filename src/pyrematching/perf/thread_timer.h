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

#ifndef PYREMATCHING_PERF_THREAD_TIMER_H
#define PYREMATCHING_PERF_THREAD_TIMER_H

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#if defined(__x86_64__) || defined(__i386__)
#include <x86intrin.h>
#define PYREMATCHING_HAVE_TSC 1
#endif

#if defined(__linux__)
#include <time.h>
#define PYREMATCHING_HAVE_THREAD_CPUTIME 1
#if defined(PYREMATCHING_HAVE_TSC)
#include <linux/perf_event.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#define PYREMATCHING_HAVE_THREAD_CYCLES 1
#endif
#endif

/// Thread-scoped stopwatches for the profile.
///
/// ## Why not `rdtsc`
///
/// The obvious "high precision" answer is `rdtsc`, and it is the wrong one for this profile. The
/// TSC is a *free-running wall clock*: it keeps counting while the thread is descheduled, so a shot
/// that loses the CPU to the scheduler is charged the scheduler's time. That is precisely the
/// instability the per-shot percentiles suffer from, and `rdtsc` does not fix it — it only makes
/// each reading cheaper. On a machine with an invariant TSC (`constant_tsc` + `nonstop_tsc`) and
/// `clocksource=tsc` there is barely even that to gain, because `steady_clock` is already served out
/// of the vDSO by the same counter: measured 22.8 ns/read against raw `rdtsc`'s 9.5 ns.
///
/// ## What is used instead
///
/// A per-thread `perf_event` cycle counter, read from user space with `rdpmc` through the event's
/// mmap page. The kernel virtualises the counter to the thread: it is saved and restored across
/// context switches, so it advances **only while this thread is on a CPU**. Measured at 8.3 ns per
/// reading on the development machine — `rdtsc`-class cost with the semantics actually wanted.
///
/// The alternative with those semantics is `clock_gettime(CLOCK_THREAD_CPUTIME_ID)`, which is a real
/// syscall (no vDSO route) at 190.8 ns per reading. That is 20x the cost, against stage timers such
/// as `dual_scan_ns` that run to tens of nanoseconds, so it is the fallback rather than the default.
///
/// ## What is being counted, exactly
///
/// **User-mode** core cycles of this thread. `exclude_kernel` is forced on — unprivileged
/// `perf_event_open` requires it at the usual `perf_event_paranoid = 2` — so time the kernel spends
/// on the thread's behalf (page faults, syscalls, the `getrusage` of `PreemptionProbe` itself) is
/// not charged to the interval. Every timer in the profile is on the same footing, so the stage
/// split and its sum stay consistent; what shifts is that `total_ns` no longer includes the
/// thread's kernel time.
///
/// Core cycles are not a fixed-rate unit — they track the core's actual frequency, so the
/// cycles-to-nanoseconds factor drifts under boost and thermal scaling. It is calibrated once,
/// lazily, on a busy loop that has first warmed the core up, and the best of three windows is kept
/// so that a descheduled calibration window cannot inflate it. Frequency drift after that point is
/// the one accuracy cost of this backend against a fixed-rate clock, and it is the reason
/// `timer_backend_name()` exists: campaigns should record which backend produced their numbers.
///
/// ## Backend selection
///
/// Chosen once per thread, on first use, most-preferred first, and overridable for A/B work with
/// the `PYREMATCHING_TIMER` environment variable (`thread` / `cputime` / `tsc` / `chrono`):
///
/// | backend       | ns/read | advances only while the thread runs |
/// |---------------|---------|-------------------------------------|
/// | `thread`      |     8.3 | yes                                 |
/// | `cputime`     |   190.8 | yes                                 |
/// | `tsc`         |     9.5 | **no**                              |
/// | `chrono`      |    22.8 | **no**                              |
///
/// The two thread-scoped backends are tried first, so a machine that refuses `perf_event_open`
/// (`perf_event_paranoid = 3`, a container without the PMU, a non-x86 host) degrades to a slower
/// thread-scoped clock rather than silently reverting to wall time. Only a non-Linux host reaches
/// the wall-clock backends by default, and `timer_backend_name()` says so.
namespace pm {
namespace perf {

enum class TimerBackend {
    /// Per-thread `perf_event` core-cycle counter, read via `rdpmc`.
    THREAD_CYCLES,
    /// `clock_gettime(CLOCK_THREAD_CPUTIME_ID)`.
    THREAD_CPUTIME,
    /// Raw `rdtsc`. Wall time: keeps counting while descheduled.
    TSC,
    /// `std::chrono::steady_clock`, named explicitly rather than through `high_resolution_clock` —
    /// the latter is a *typedef* for `system_clock` on libstdc++, which is not monotonic and can
    /// step backwards under NTP. Wall time.
    STEADY_CLOCK,
};

inline const char* backend_name(TimerBackend backend) {
    switch (backend) {
        case TimerBackend::THREAD_CYCLES:
            return "thread_cycles_rdpmc";
        case TimerBackend::THREAD_CPUTIME:
            return "clock_thread_cputime_id";
        case TimerBackend::TSC:
            return "rdtsc";
        case TimerBackend::STEADY_CLOCK:
            return "steady_clock";
    }
    return "unknown";
}

/// True iff the backend stops while the thread is off the CPU. Reported beside timings so that a
/// campaign which fell back to a wall-clock backend cannot be read as if it had not.
inline bool backend_is_thread_scoped(TimerBackend backend) {
    return backend == TimerBackend::THREAD_CYCLES || backend == TimerBackend::THREAD_CPUTIME;
}

namespace internal {

/// Deliberately a POD with no destructor. It owns a file descriptor and a mapped page which are
/// **not** released at thread exit: the timers are used from static-lifetime objects and from
/// destructors, and tearing the clock down underneath them would be worse than holding one fd and
/// one page for the life of the (single-threaded, by design) process.
struct ThreadClock {
    bool initialized;
    TimerBackend backend;
    double ns_per_tick;
#if defined(PYREMATCHING_HAVE_THREAD_CYCLES)
    struct perf_event_mmap_page* page;
#endif
};

inline void compiler_barrier() {
    asm volatile("" ::: "memory");
}

#if defined(PYREMATCHING_HAVE_THREAD_CYCLES)
/// The seqlock read of the event's mmap page, as `perf_event.h` documents it. `lock` is bumped by
/// the kernel around every schedule-in/out of the event, so a reading that straddles one is retried.
inline uint64_t read_perf_cycles(struct perf_event_mmap_page* page) {
    uint64_t count;
    uint32_t seq;
    do {
        seq = page->lock;
        compiler_barrier();
        uint32_t index = page->index;
        count = page->offset;
        if (index == 0) {
            // The event is not currently scheduled onto a PMU counter, so there is nothing to add:
            // `offset` already holds everything counted up to the last schedule-out. For a pinned
            // self-monitoring event read by the very thread it counts this does not arise, but a
            // stale-but-monotone reading beats a garbage one if it ever does.
            break;
        }
        int64_t pmc = (int64_t)__rdpmc(index - 1);
        // Sign-extend from the counter's width so that a wrap since `offset` was written subtracts
        // rather than adding a spurious 2^width.
        pmc <<= 64 - page->pmc_width;
        pmc >>= 64 - page->pmc_width;
        count += (uint64_t)pmc;
        compiler_barrier();
    } while (page->lock != seq);
    return count;
}
#endif

inline uint64_t steady_clock_ns() {
    return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
               std::chrono::steady_clock::now().time_since_epoch())
        .count();
}

#if defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
inline uint64_t thread_cputime_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
}
#endif

/// Ticks are cycles on `THREAD_CYCLES` and on `TSC`, and already nanoseconds on the other two.
/// `ThreadClock::ns_per_tick` reconciles them.
inline uint64_t read_ticks(const ThreadClock& clock) {
    switch (clock.backend) {
#if defined(PYREMATCHING_HAVE_THREAD_CYCLES)
        case TimerBackend::THREAD_CYCLES:
            return read_perf_cycles(clock.page);
#endif
#if defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
        case TimerBackend::THREAD_CPUTIME:
            return thread_cputime_ns();
#endif
#if defined(PYREMATCHING_HAVE_TSC)
        case TimerBackend::TSC:
            return __rdtsc();
#endif
        default:
            return steady_clock_ns();
    }
}

/// Nanoseconds per tick, for the backends whose ticks are not already nanoseconds.
///
/// The reference clock must have the **same scope as the counter being calibrated**, which is the
/// whole trick here. Calibrating the thread-scoped cycle counter against wall time would be a race
/// between the two: a calibration window the thread is descheduled during accrues wall nanoseconds
/// without accruing thread cycles, silently inflating ns-per-cycle and stretching every later
/// measurement. Against `CLOCK_THREAD_CPUTIME_ID` both sides stop together, so the ratio is
/// unaffected by however much the machine interferes — its 190 ns cost does not matter for four
/// readings taken once. The TSC backend is wall-clock, so it calibrates against `steady_clock`.
///
/// Busy-waits rather than sleeps: sleeping would measure the scheduler, and on a thread-scoped
/// counter it would not advance the counter at all. The first window is discarded as a 25 ms
/// warm-up, because a core that starts idle needs tens of milliseconds to reach its boost clock and
/// a shorter ramp measured `ns_per_tick` as much as 8% high on this machine — the counter is
/// cycles, so calibrating at an unboosted frequency stretches every subsequent reading.
///
/// The **smallest** of the kept windows wins, for the same reason: a window that caught a
/// frequency excursion can only have counted fewer cycles per nanosecond, never more, so the
/// minimum is the least-contaminated estimate. It cannot be gamed by descheduling the way a
/// wall-clock reference could, because the reference here stops with the counter.
///
/// What survives is a systematic bias, not a variance: calibration runs one thread, so it sees
/// single-core boost, while a long campaign settles at the lower all-core sustained clock. That
/// makes absolute nanoseconds read a little low under sustained load. Ratios between two timers in
/// the same run — which is what the speedup tables are built from — divide it out entirely.
inline double calibrate_ns_per_tick(const ThreadClock& clock, uint64_t (*reference_ns)()) {
    double best = 0;
    for (int window = 0; window < 4; window++) {
        uint64_t wall_start = steady_clock_ns();
        uint64_t ref_start = reference_ns();
        uint64_t tick_start = read_ticks(clock);
        uint64_t target = (window == 0 ? 25u : 10u) * 1000000ull;
        while (steady_clock_ns() - wall_start < target) {
        }
        uint64_t ticks = read_ticks(clock) - tick_start;
        uint64_t ns = reference_ns() - ref_start;
        if (window == 0 || ticks == 0)
            continue;  // warm-up, or a window that counted nothing
        double ratio = (double)ns / (double)ticks;
        if (best == 0 || ratio < best)
            best = ratio;
    }
    return best;
}

#if defined(PYREMATCHING_HAVE_THREAD_CYCLES)
/// Opens a pinned, user-mode, per-thread cycle counter and maps its control page. Returns null if
/// the kernel refuses the event, if there is no PMU to put it on, or if this kernel will not let
/// user space read the counter with `rdpmc`.
inline struct perf_event_mmap_page* open_thread_cycle_counter() {
    struct perf_event_attr attr;
    memset(&attr, 0, sizeof(attr));
    attr.type = PERF_TYPE_HARDWARE;
    attr.size = sizeof(attr);
    attr.config = PERF_COUNT_HW_CPU_CYCLES;
    // Required at the usual `perf_event_paranoid = 2`; see the header comment on what this excludes.
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    // Never multiplex this event. A multiplexed counter reports a *sample* scaled by
    // time_enabled/time_running, which would silently undercount every interval here.
    attr.pinned = 1;
    // pid 0 = the calling thread, cpu -1 = wherever it runs. This is what makes the counter follow
    // the thread across CPUs and stop while it is off them.
    long fd = syscall(__NR_perf_event_open, &attr, /*pid=*/0, /*cpu=*/-1, /*group_fd=*/-1, /*flags=*/0);
    if (fd < 0)
        return nullptr;
    void* mapped = mmap(nullptr, (size_t)sysconf(_SC_PAGESIZE), PROT_READ, MAP_SHARED, (int)fd, 0);
    if (mapped == MAP_FAILED) {
        close((int)fd);
        return nullptr;
    }
    auto* page = (struct perf_event_mmap_page*)mapped;
    if (!page->cap_user_rdpmc || page->index == 0) {
        munmap(mapped, (size_t)sysconf(_SC_PAGESIZE));
        close((int)fd);
        return nullptr;
    }
    // The fd is deliberately not closed: closing it destroys the event and the mapping goes dead.
    return page;
}
#endif

/// `PYREMATCHING_TIMER`, or null. Forces a backend. An unrecognised value is ignored rather than
/// honoured: a typo must not be able to silently downgrade a campaign to a wall-clock timer, which
/// is the one failure here that produces plausible numbers instead of obviously broken ones.
inline const char* backend_override() {
    const char* value = getenv("PYREMATCHING_TIMER");
    if (value == nullptr || value[0] == '\0')
        return nullptr;
    if (strcmp(value, "thread") == 0 || strcmp(value, "cputime") == 0 || strcmp(value, "tsc") == 0 ||
        strcmp(value, "chrono") == 0) {
        return value;
    }
    return nullptr;
}

inline void initialize(ThreadClock& clock) {
    const char* forced = backend_override();
    bool want = forced == nullptr;

#if defined(PYREMATCHING_HAVE_THREAD_CYCLES)
    if (want || strcmp(forced, "thread") == 0) {
        clock.page = open_thread_cycle_counter();
        if (clock.page != nullptr) {
            clock.backend = TimerBackend::THREAD_CYCLES;
            clock.ns_per_tick = calibrate_ns_per_tick(clock, thread_cputime_ns);
            if (clock.ns_per_tick > 0) {
                clock.initialized = true;
                return;
            }
        }
    }
#endif
#if defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
    if (want || strcmp(forced, "cputime") == 0) {
        clock.backend = TimerBackend::THREAD_CPUTIME;
        clock.ns_per_tick = 1.0;
        clock.initialized = true;
        return;
    }
#endif
#if defined(PYREMATCHING_HAVE_TSC)
    if (forced != nullptr && strcmp(forced, "tsc") == 0) {
        clock.backend = TimerBackend::TSC;
        clock.ns_per_tick = calibrate_ns_per_tick(clock, steady_clock_ns);
        if (clock.ns_per_tick > 0) {
            clock.initialized = true;
            return;
        }
    }
#endif
    clock.backend = TimerBackend::STEADY_CLOCK;
    clock.ns_per_tick = 1.0;
    clock.initialized = true;
}

inline ThreadClock& thread_clock() {
    // Zero-initialised, trivially destructible, so this is a TLS address plus a predictable branch
    // on the hot path rather than a guarded function-local static with a destructor to register.
    static thread_local ThreadClock clock{};
    if (!clock.initialized)
        initialize(clock);
    return clock;
}

}  // namespace internal

/// Which backend this thread settled on. Record it beside any timing a campaign publishes.
inline TimerBackend current_backend() {
    return internal::thread_clock().backend;
}

inline const char* current_backend_name() {
    return backend_name(current_backend());
}

/// A stopwatch over thread run time. See the header comment for what "run time" excludes.
///
/// Not copied between threads: the tick origin is this thread's counter.
struct ThreadTimer {
    uint64_t start_tick{0};

    inline void start() {
        start_tick = internal::read_ticks(internal::thread_clock());
    }

    inline long long elapsed_ns() const {
        const internal::ThreadClock& clock = internal::thread_clock();
        uint64_t ticks = internal::read_ticks(clock) - start_tick;
        return (long long)((double)ticks * clock.ns_per_tick);
    }
};

}  // namespace perf
}  // namespace pm

#endif  // PYREMATCHING_PERF_THREAD_TIMER_H
