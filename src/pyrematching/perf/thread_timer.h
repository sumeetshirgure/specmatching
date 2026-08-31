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
#include <cstdio>
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
#elif defined(__APPLE__)
#include <time.h>
// `CLOCK_THREAD_CPUTIME_ID` arrived with the 10.12 SDK, as a macro over an enumerator, so testing the
// macro is testing whether this SDK has the clock at all. Darwin serves it from the
// `thread_selfusage` mach trap — this thread's user + kernel time, in nanoseconds, saved and restored
// across context switches — which is the same semantics Linux's `CLOCK_THREAD_CPUTIME_ID` has and
// exactly what the profile wants. Whether the *runtime* honours it is settled by
// `thread_cputime_available()`, not here: an SDK that declares it is not a kernel that serves it.
#if defined(CLOCK_THREAD_CPUTIME_ID)
#define PYREMATCHING_HAVE_THREAD_CPUTIME 1
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
/// | backend       | ns/read | advances only while the thread runs | where                 |
/// |---------------|---------|-------------------------------------|-----------------------|
/// | `thread`      |     8.3 | yes                                 | Linux + x86 PMU       |
/// | `cputime`     |   190.8 | yes                                 | Linux, macOS 10.12+   |
/// | `tsc`         |     9.5 | **no**                              | x86, forced only      |
/// | `chrono`      |    22.8 | **no**                              | anywhere              |
///
/// The two thread-scoped backends are tried first, so a machine that refuses `perf_event_open`
/// (`perf_event_paranoid = 3`, a container without the PMU, a non-x86 host) degrades to a slower
/// thread-scoped clock rather than silently reverting to wall time.
///
/// ## macOS
///
/// There is no `perf_event_open` there and no unprivileged route to the PMU, so `thread` is not
/// available and `cputime` is the backend a macOS run should land on:
/// `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` is served by the `thread_selfusage` mach trap, counts
/// this thread's user + kernel nanoseconds, and stops while the thread is off a CPU — the same
/// contract the Linux fallback has. It is *not* `mach_absolute_time`, which is a wall clock and
/// would have exactly the defect this whole header exists to avoid. Two consequences worth knowing
/// when reading a macOS campaign against a Linux one:
///
///   - the reading is a trap rather than a vDSO call, so it is the ~200 ns-class backend, not the
///     8 ns one. Stage timers in the tens of nanoseconds are dominated by it; sums over a shot are
///     not. `timer_backend_name()` is in every log header for precisely this comparison;
///   - `cputime` charges kernel time to the interval and `thread` does not (it is forced to
///     `exclude_kernel`), so a macOS `total_ns` includes syscall and fault time that the same run on
///     Linux would leave out.
///
/// A macOS build only reaches the wall-clock backends if the runtime refuses the clock outright, and
/// `backend_is_thread_scoped()` — reported in every log header, and warned about on stderr by
/// `warn_if_wall_clock` — is what says so rather than leaving it to be assumed.
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

/// Does the *runtime* actually serve the clock the SDK declared?
///
/// Asked once, and asked at all because of macOS: the clock is a compile-time feature of the SDK and
/// a run-time feature of the kernel, and a binary built on 10.12+ can be run somewhere the trap is
/// refused. Two readings rather than one, because the failure that matters is not only `EINVAL` —
/// a clock that returns success and never advances is a wall-clock fallback wearing a thread-scoped
/// name, and a busy loop between the two readings is the cheapest way to catch it. The loop is
/// bounded by iterations, not by time, so this cannot hang on a stopped clock.
inline bool thread_cputime_available() {
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) != 0)
        return false;
    uint64_t before = (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
    volatile uint64_t sink = 0;
    // Rounds of ~0.5 ms of work until the clock moves, up to ~30 ms of CPU in all. Bounded by
    // iterations rather than by a deadline so a stopped clock ends the loop instead of owning it,
    // and given that many rounds so that a coarse-but-working clock is not mistaken for a stopped
    // one — a clock that cannot move in 30 ms of CPU time is no use for timing microseconds anyway.
    for (int round = 0; round < 64; round++) {
        for (int i = 0; i < 200000; i++)
            sink = sink + (uint64_t)i;
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) != 0)
            return false;
        uint64_t after = (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
        if (after > before)
            return true;
    }
    return false;
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
/// The busy-wait must spend the window in **user mode**, which is why it is a batch of arithmetic
/// between clock reads rather than a spin on the clock itself. `THREAD_CYCLES` counts user cycles
/// only (`exclude_kernel` is forced), while every reference clock counts kernel time too, so any
/// kernel time inside the window lands in the numerator and not in the denominator and inflates
/// `ns_per_tick` directly. A tight `while (steady_clock_ns() - start < target) {}` is the worst case
/// of exactly that: where the clocksource is not vDSO-able — `hpet`, `acpi_pm`, a VM without a
/// paravirtual clock — every iteration is a syscall, the window is ~95% kernel time, and the
/// calibration comes out **60x** high. That was measured here, on an `hpet` box where a
/// `steady_clock` read costs 1.3 us: 15.3 ns/cycle against a true 0.253, i.e. every latency in every
/// report inflated 60-fold while the ratios between them stayed entirely plausible. Batching the
/// window into millisecond runs of arithmetic between clock reads brought the same host to 0.2575,
/// 1.6% high — which is the frequency bias the paragraph below is about, and no longer a broken
/// instrument. `calibration_looks_sane` is the backstop for whatever this reasoning has still
/// missed.
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
    // Volatile so the batch below is real work rather than something the optimiser folds away, and
    // hoisted out of the loop so the whole calibration touches one cache line.
    volatile uint64_t sink = 0;
    for (int window = 0; window < 4; window++) {
        uint64_t wall_start = steady_clock_ns();
        uint64_t ref_start = reference_ns();
        uint64_t tick_start = read_ticks(clock);
        uint64_t target = (window == 0 ? 25u : 10u) * 1000000ull;
        do {
            // Sub-millisecond of user-mode arithmetic per clock read, so the window is user mode to
            // within a fraction of a percent even where a `steady_clock` read costs the 1.3 us
            // measured on the `hpet` host above. Overshooting the target by up to one batch is
            // harmless: the window is measured, not assumed.
            for (int i = 0; i < 1000000; i++)
                sink = sink + (uint64_t)i;
        } while (steady_clock_ns() - wall_start < target);
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

/// Is a calibrated cycle counter's ns-per-cycle a number a CPU could actually have?
///
/// The backstop on everything `calibrate_ns_per_tick` reasons about. A miscalibrated cycle counter
/// does not fail, it *scales*: every stage, every mean and every tail moves by the same factor, the
/// ratios between them stay exactly as plausible as before, and nothing in a report looks wrong. The
/// 60x seen on an `hpet` host was found by comparing backends, not by reading a number that looked
/// odd. So the range is checked rather than trusted, and a counter outside it is refused — falling
/// back to a slower thread-scoped clock that is right beats keeping a fast one that is wrong.
///
/// The bounds are deliberately loose: `[0.02, 2.0]` ns/cycle is 0.5 GHz to 50 GHz, which no real core
/// leaves and no plausible future one will either. This is a sanity check on the *instrument*, not a
/// judgement about the machine.
inline bool calibration_looks_sane(double ns_per_tick) {
    return ns_per_tick >= 0.02 && ns_per_tick <= 2.0;
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
            // A cycle counter is only as good as its calibration, and a bad one is silent — see
            // `calibration_looks_sane`. Refused rather than kept, so the run lands on `cputime` and
            // says `clock_thread_cputime_id` in its log header.
            if (calibration_looks_sane(clock.ns_per_tick)) {
                clock.initialized = true;
                return;
            }
        }
    }
#endif
#if defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
    // The only thread-scoped backend macOS has, and the one a macOS run is expected to land on.
    // Probed rather than assumed: on Linux the syscall is always there, on macOS the SDK declaring
    // the clock is not the kernel serving it, and falling through to a wall clock while still
    // *claiming* `clock_thread_cputime_id` would be the one failure mode that produces plausible
    // numbers instead of obviously broken ones.
    if (want || strcmp(forced, "cputime") == 0) {
        if (thread_cputime_available()) {
            clock.backend = TimerBackend::THREAD_CPUTIME;
            clock.ns_per_tick = 1.0;
            clock.initialized = true;
            return;
        }
    }
#endif
#if defined(PYREMATCHING_HAVE_TSC)
    if (forced != nullptr && strcmp(forced, "tsc") == 0) {
        clock.backend = TimerBackend::TSC;
        clock.ns_per_tick = calibrate_ns_per_tick(clock, steady_clock_ns);
        // The TSC is a fixed-rate counter rather than a core-cycle one, so this is a wider net than
        // it looks — but a TSC calibrated against a wall clock is measuring like against like, and
        // the check is here for the same reason as above: a scale error would be invisible.
        if (calibration_looks_sane(clock.ns_per_tick)) {
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

/// Wall time, for the one job that needs a wall clock: measuring how much of an interval this thread
/// spent *not* running (`PreemptionProbe`). Never charge a stage with it.
inline uint64_t wall_clock_ns() {
    return internal::steady_clock_ns();
}

/// This thread's CPU time — user + kernel — or false where the platform has no thread-scoped clock.
///
/// Deliberately not routed through `ThreadTimer`: this is the *reference* clock, the one whose gap
/// to wall time says the thread was descheduled, so it must be the thread-scoped clock itself and
/// not whichever backend the timer settled on.
inline bool thread_cpu_ns(uint64_t& out) {
#if defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
    struct timespec ts;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) != 0)
        return false;
    out = (uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec;
    return true;
#else
    (void)out;
    return false;
#endif
}

/// Says, loudly and on `stream`, that this run's numbers are wall-clock numbers. Returns whether it
/// warned, so a caller can record the fact as well as print it.
///
/// Worth a banner rather than a parenthesis in a status line because of what a wall-clock fallback
/// does to a campaign: it does not fail, it does not look wrong, it just quietly charges every stage
/// for whatever else the machine was doing, and every mean and every speedup built on it inherits
/// that. The failure is silent by nature, so the warning cannot be.
inline bool warn_if_wall_clock(std::FILE* stream) {
    if (backend_is_thread_scoped(current_backend()))
        return false;
    std::fprintf(
        stream,
        "\n"
        "  ##########################################################################\n"
        "  WARNING: latency is being measured with a WALL CLOCK, backend '%s'.\n"
        "  ##########################################################################\n"
        "  No thread-scoped clock was available, so every timing on this run includes\n"
        "  whatever time the thread spent off the CPU: means, tails and speedups all\n"
        "  carry the scheduler's interference inside them and are upper bounds only.\n"
        "  Linux: needs perf_event_open (perf_event_paranoid <= 2) or clock_gettime\n"
        "         (CLOCK_THREAD_CPUTIME_ID).\n"
        "  macOS: needs clock_gettime(CLOCK_THREAD_CPUTIME_ID), i.e. macOS 10.12+.\n"
        "  Check PYREMATCHING_TIMER, which forces the backend and may be set to chrono.\n"
        "\n",
        current_backend_name());
    std::fflush(stream);
    return true;
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
