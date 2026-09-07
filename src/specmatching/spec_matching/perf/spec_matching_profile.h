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

#ifndef SPECMATCHING_SPEC_MATCHING_PERF_SPEC_MATCHING_PROFILE_H
#define SPECMATCHING_SPEC_MATCHING_PERF_SPEC_MATCHING_PROFILE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "specmatching/perf/thread_timer.h"
#include "specmatching/sparse_blossom/ints.h"
#include "specmatching/spec_matching/truncation/harvest.h"

#if defined(__linux__)
#include <sys/resource.h>
#endif

/// Linux hands out a per-thread context-switch counter and nothing else does, so `RUSAGE_THREAD` is
/// the preferred mechanism where it exists and `SPECMATCHING_PREEMPTION_CLOCK_GAP` is the fallback
/// used everywhere else that has a thread-scoped clock — macOS in particular. See `PreemptionProbe`.
///
/// `-DSPECMATCHING_FORCE_CLOCK_GAP_PROBE` selects the fallback on a machine that has the counters,
/// which is how the macOS path is exercised from a Linux CI box — the mechanism is chosen at compile
/// time, so without it there is no way to run that code at all.
#if defined(SPECMATCHING_FORCE_CLOCK_GAP_PROBE) && defined(SPECMATCHING_HAVE_THREAD_CPUTIME)
#define SPECMATCHING_PREEMPTION_CLOCK_GAP 1
#elif defined(__linux__) && defined(RUSAGE_THREAD)
#define SPECMATCHING_PREEMPTION_SWITCH_COUNTERS 1
#elif defined(SPECMATCHING_HAVE_THREAD_CPUTIME)
#define SPECMATCHING_PREEMPTION_CLOCK_GAP 1
#endif

namespace pm {
namespace spec_matching {

/// The stopwatch every stage of the profile is timed with (§M6.1 timer backends, §M6.4 timer
/// discipline). Thread-scoped: it advances only while this thread is on a CPU, so a shot the
/// scheduler interrupts is no longer charged for the time it spent off the CPU. See
/// `perf/thread_timer.h` for the backends, their measured per-reading costs, and what "thread run
/// time" excludes.
using HiResTimer = pm::perf::ThreadTimer;

/// Was this thread descheduled while the shot was being timed? (§M6.4 timer discipline.)
///
/// `HiResTimer` no longer *charges* the shot for time spent off the CPU, so this is no longer what
/// stands between a percentile and the scheduler. It is kept because being descheduled costs the
/// shot more than the wall time it loses — the caches and branch predictors it resumes on are not
/// the ones it left — so an interrupted shot is still an outlier, just a far smaller one. It is
/// also the cross-check on the clock: on a thread-scoped backend a contaminated shot should now
/// look much like its neighbours, and if it does not, the backend is not doing what it claims.
///
/// ## How it is detected, and why that differs by platform
///
/// `getrusage(RUSAGE_THREAD)` counts this thread's voluntary and involuntary context switches, so on
/// Linux the question is answered exactly: the counters moved, or they did not.
///
/// macOS has no `RUSAGE_THREAD` and exposes no per-thread switch counter at all — not through
/// `getrusage`, not through `proc_pid_rusage`, not through `thread_info` — so the same question is
/// answered by **the gap between two clocks** instead. Over any interval, wall time minus this
/// thread's CPU time is the time the thread was not running; a run of that anywhere near the cost of
/// a context switch means the thread lost the CPU. This is a different instrument from the Linux one
/// and it is worth being clear about how it differs:
///
///   - it measures *how much* time was lost rather than *how many* switches happened, so a switch
///     that returns the CPU inside the slack below is not seen. That is the intended trade: the
///     reason a contaminated shot is dropped is the state it lost, and a switch too short to show up
///     in the clock gap did not have time to lose much of it;
///   - it charges the same way for anything else that takes the thread off the CPU — a page fault
///     that goes to disk, a blocking syscall — which for this purpose is a feature, since those cost
///     the shot the same caches;
///   - both readings are taken inside `sample()`, in a fixed order, so the skew between them cancels
///     between two samples to within the variation of two clock reads (tens of nanoseconds against a
///     two-microsecond slack).
///
/// The slack is `SPECMATCHING_PREEMPTION_SLACK_NS`, default 2000. Well above the noise of four clock
/// readings, well below the tens of microseconds a real deschedule costs.
///
/// Run against the switch counters on the same intervals on a Linux box (which is what
/// `SPECMATCHING_FORCE_CLOCK_GAP_PROBE` is for), the gap flags a **superset**: over 200 intervals of
/// ~30 us it flagged 35 where the counters flagged 29, and every one of those 29 was among them. The
/// extras are intervals that lost 2 us or more without a thread context switch being charged for it,
/// interrupt handling most of them — which is off-CPU time that costs the shot its caches just the
/// same. Over intervals of a few microseconds, the length a shot actually is, the mean gap on an idle
/// machine measured 0.8 us against the 2 us slack, so the skew between the two clocks does not flag
/// anything on its own.
///
/// Where neither mechanism exists the probe reports "not contaminated" and `contaminated_shot_rate`
/// reads 0, which the artifact says explicitly rather than implying the machine was quiet. Read
/// `mechanism_name()` beside the rate, for the same reason the timer backend is read beside the
/// timings: a rate of 0 from `none` is not a measurement.
struct PreemptionProbe {
    /// Switch counters, on the platforms that have them.
    long voluntary{0};
    long involuntary{0};
    /// The two clocks whose divergence stands in for them where they do not.
    uint64_t wall_ns{0};
    uint64_t cpu_ns{0};
    bool supported{false};

    inline void sample() {
#if defined(SPECMATCHING_PREEMPTION_SWITCH_COUNTERS)
        struct rusage usage;
        if (getrusage(RUSAGE_THREAD, &usage) == 0) {
            voluntary = usage.ru_nvcsw;
            involuntary = usage.ru_nivcsw;
            supported = true;
            return;
        }
#elif defined(SPECMATCHING_PREEMPTION_CLOCK_GAP)
        // Order matters and is fixed: CPU first, wall second, in both samples. The interval's wall
        // time then includes the cost of one CPU-clock read at each end and the interval's CPU time
        // does not, which biases the gap by a constant that is the same on every shot and is an order
        // of magnitude under the slack.
        uint64_t cpu = 0;
        if (pm::perf::thread_cpu_ns(cpu)) {
            cpu_ns = cpu;
            wall_ns = pm::perf::wall_clock_ns();
            supported = true;
            return;
        }
#endif
        supported = false;
    }

    /// True when this thread lost the CPU between `before` and this sample.
    inline bool switched_since(const PreemptionProbe& before) const {
        if (!supported || !before.supported)
            return false;
#if defined(SPECMATCHING_PREEMPTION_SWITCH_COUNTERS)
        return voluntary != before.voluntary || involuntary != before.involuntary;
#elif defined(SPECMATCHING_PREEMPTION_CLOCK_GAP)
        // Unsigned throughout, so a wall clock that failed to advance while the CPU clock did cannot
        // wrap into a huge positive gap and flag every shot.
        if (wall_ns <= before.wall_ns || cpu_ns < before.cpu_ns)
            return false;
        uint64_t wall = wall_ns - before.wall_ns;
        uint64_t cpu = cpu_ns - before.cpu_ns;
        return wall > cpu && wall - cpu > slack_ns();
#else
        return false;
#endif
    }

    /// Which of the two instruments above filled this build's samples. Recorded beside the rate.
    static inline const char* mechanism_name() {
#if defined(SPECMATCHING_PREEMPTION_SWITCH_COUNTERS)
        return "rusage_thread_switches";
#elif defined(SPECMATCHING_PREEMPTION_CLOCK_GAP)
        return "wall_minus_thread_cpu";
#else
        return "none";
#endif
    }

    /// Nanoseconds of off-CPU time an interval is allowed before the shot is called contaminated.
    /// Meaningless under the switch-counter mechanism, which does not measure a duration.
    static inline uint64_t slack_ns() {
        static const uint64_t slack = [] {
            const char* value = std::getenv("SPECMATCHING_PREEMPTION_SLACK_NS");
            if (value == nullptr || value[0] == '\0')
                return (uint64_t)2000;
            char* end = nullptr;
            unsigned long long parsed = std::strtoull(value, &end, 10);
            // A typo must not silently widen the slack until nothing is ever flagged.
            if (end == value || *end != '\0')
                return (uint64_t)2000;
            return (uint64_t)parsed;
        }();
        return slack;
    }
};

/// §3.5.2's distributions, over one shot's components, `H` edges and defects.
///
/// Kept apart from the per-shot profile because their unit is the component, the edge or the node
/// rather than the shot: the profile stays a flat row of scalars the pybind layer can column-ise,
/// and these are folded straight into the campaign accumulator instead. One instance is reused
/// across shots, so a profiled shot allocates nothing for them either.
///
/// **Component size is the only per-component distribution.** The degree histogram and the two
/// `H`-subgraph diameter histograms are gone, along with the adjacency and the all-pairs walks that
/// filled them.
///
/// **Unit.** The two weight histograms are binned in *sixteenths of the horizon* — bin `k` counts
/// values in `[k * T / 16, (k + 1) * T / 16)`, and the last bin overflows — so both share one axis
/// and a plot can label it without knowing the normalising constant. `H`'s defect-defect edges are
/// `<= 2T` by construction, which is exactly the last bin of `edge_weight_hist`.
///
/// `edge_weight_hist` and `bcost_hist` are the two §C.2 accumulators `sparse_graph_stats`
/// deliberately does not write (§3.5.2, "Not collected"). They survived the cut to size-only
/// component statistics because neither is keyed by component — one is per `H` edge, the other per
/// defect — and the latency profiler still writes them.
///
/// Filled outside every timed region, and never counted toward any reported latency.
struct ComponentHistograms {
    static constexpr size_t BINS_PER_T = 16;
    /// Component sizes `0..31`, last bin overflowing. Bin 0 is always empty: a component has at
    /// least one member. `configure` moves this; `configure_uncapped` removes it.
    static constexpr size_t SIZE_HIST_BINS = 33;
    /// `w_int` over `H`'s defect-defect edges: `0 .. 2T`, which is the whole range.
    static constexpr size_t WEIGHT_HIST_BINS = 33;
    /// `bcost_int` over defects that have a boundary entry within `R`. Values past `2T` — legal,
    /// since `R >= 2 * T_max` — land in the overflow bin.
    static constexpr size_t BCOST_HIST_BINS = 33;

    std::vector<uint64_t> size_hist = std::vector<uint64_t>(SIZE_HIST_BINS, 0);
    std::vector<uint64_t> edge_weight_hist = std::vector<uint64_t>(WEIGHT_HIST_BINS, 0);
    std::vector<uint64_t> bcost_hist = std::vector<uint64_t>(BCOST_HIST_BINS, 0);

    /// Whether `size_hist` has an overflow bin at all.
    ///
    /// Set by `configure_uncapped`, and false everywhere else. It governs `size_hist` alone; the two
    /// weight-binned histograms are ranged in multiples of the horizon, which is not a
    /// component-size cap, and their last bin still overflows either way.
    bool uncapped{false};

    /// Moves the component-size cap. The cap is the **overflow bin**, so a cap of `c` gives `c + 1`
    /// bins and everything at or above `c` lands in the last one — the `residual_size_hist`
    /// convention. Resizes and zeroes; two histogram sets that are added together must have been
    /// configured the same way.
    void configure(size_t size_cap) {
        uncapped = false;
        size_hist.assign(size_cap + 1, 0);
    }

    /// Removes the component-size cap: `size_hist` grows to fit whatever it is handed, so bin `k`
    /// counts exactly the size `k` at every `k` and there is no overflow bin to read past.
    ///
    /// This is what `sparse_graph_stats` configures. That binary's entire output is the component
    /// size distribution, and a size above which components silently pile into one terminal bin is a
    /// distribution with its tail cut off — the part the escalation rate is actually about.
    /// Everything on a decode path keeps `configure`'s fixed bins, where the array is sized once and
    /// never reallocates.
    void configure_uncapped() {
        uncapped = true;
        // One bin, not none: an empty array would make `add` on an all-empty cell a special case
        // for no gain.
        size_hist.assign(1, 0);
    }

    /// Counts one component into `size_hist`: clamped into the overflow bin when capped, growing the
    /// array when not.
    inline void add_component_size(uint64_t size) {
        if (uncapped) {
            if (size >= size_hist.size())
                size_hist.resize((size_t)size + 1, 0);
            size_hist[(size_t)size]++;
            return;
        }
        size_hist[count_bin(size, size_hist.size())]++;
    }

    /// Zeroes every bin without touching the buffers, so the next shot reuses the same storage.
    void clear() {
        std::fill(size_hist.begin(), size_hist.end(), (uint64_t)0);
        std::fill(edge_weight_hist.begin(), edge_weight_hist.end(), (uint64_t)0);
        std::fill(bcost_hist.begin(), bcost_hist.end(), (uint64_t)0);
    }

    /// The bin a weight-like quantity falls in, in sixteenths of `horizon`. A non-positive horizon
    /// cannot happen on a decoding path — it is rejected at construction — and reads as bin 0 rather
    /// than dividing by it.
    static inline size_t weight_bin(pm::cumulative_time_int value, pm::cumulative_time_int horizon, size_t bins) {
        if (horizon <= 0 || value <= 0)
            return 0;
        pm::cumulative_time_int bin = value * (pm::cumulative_time_int)BINS_PER_T / horizon;
        return (size_t)std::min<pm::cumulative_time_int>(bin, (pm::cumulative_time_int)bins - 1);
    }

    static inline size_t count_bin(uint64_t value, size_t bins) {
        return (size_t)std::min<uint64_t>(value, (uint64_t)bins - 1);
    }

    void add(const ComponentHistograms& other) {
        auto add_bins = [](std::vector<uint64_t>& into, const std::vector<uint64_t>& from) {
            // Two sets configured differently would silently drop or misalign bins, and the range is
            // a run-level choice that cannot change mid-campaign.
            assert(into.size() == from.size() && "histogram sets were configured with different ranges");
            for (size_t i = 0; i < into.size(); i++)
                into[i] += from[i];
        };
        // Uncapped, two size histograms may legitimately differ in length — the one that saw the
        // larger component has the longer array, and the shorter one is that array's prefix, because
        // bin `k` means the size `k` in both. Capped, they still have to agree bin for bin.
        assert(uncapped == other.uncapped && "an uncapped histogram set was added to a capped one");
        assert(
            (uncapped || size_hist.size() == other.size_hist.size()) &&
            "histogram sets were configured with different caps");
        if (other.size_hist.size() > size_hist.size())
            size_hist.resize(other.size_hist.size(), 0);
        for (size_t i = 0; i < other.size_hist.size(); i++)
            size_hist[i] += other.size_hist[i];
        add_bins(edge_weight_hist, other.edge_weight_hist);
        add_bins(bcost_hist, other.bcost_hist);
    }
};

/// §2.6/§3.5.2's joint table: how a component's `COMPLETE`/`TRUNCATED` status is distributed over
/// its size.
///
/// This is the table the whole branch exists to make readable — the status is now a per-component
/// fact rather than a per-shot one, because every component is decided by its own truncated solve
/// (§2). Size is the only key it is ever built on; the hop-diameter twin went with the diameters.
struct ComponentStatusTable {
    /// Column 0 is `COMPLETE`, column 1 is `TRUNCATED`.
    static constexpr size_t STATUSES = 2;
    static constexpr size_t COMPLETE = 0;
    static constexpr size_t TRUNCATED = 1;

    /// Row-major `[bin][status]`, with the last bin overflowing unless `uncapped`.
    std::vector<uint64_t> counts;
    size_t bins{0};
    /// Whether the table grows to fit its key instead of clamping into a last bin. See
    /// `ComponentHistograms::uncapped`, which `sparse_graph_stats` sets alongside this one — the
    /// table is keyed by component size, so an overflow bin here cuts the tail off the same
    /// distribution the size histogram reports.
    bool uncapped{false};

    void configure(size_t bin_count) {
        uncapped = false;
        bins = bin_count;
        counts.assign(bin_count * STATUSES, 0);
    }
    static ComponentStatusTable with_bins(size_t bin_count) {
        ComponentStatusTable table;
        table.configure(bin_count);
        return table;
    }
    /// Removes the cap: bin `k` is exactly the size `k`, and rows appear as sizes are seen.
    void configure_uncapped() {
        configure(1);
        uncapped = true;
    }
    static ComponentStatusTable growing() {
        ComponentStatusTable table;
        table.configure_uncapped();
        return table;
    }
    inline void add(size_t bin, size_t status) {
        if (uncapped) {
            // Rows are appended and the array is row-major, so a resize leaves every existing
            // `[bin][status]` at the index it already had.
            if (bin >= bins) {
                bins = bin + 1;
                counts.resize(bins * STATUSES, 0);
            }
        } else {
            if (bins == 0)
                return;
            bin = std::min(bin, bins - 1);
        }
        counts[bin * STATUSES + status]++;
    }
    inline uint64_t at(size_t bin, size_t status) const {
        return counts[bin * STATUSES + status];
    }
    void add(const ComponentStatusTable& other) {
        assert(uncapped == other.uncapped && "an uncapped status table was added to a capped one");
        assert((uncapped || bins == other.bins) && "status tables were configured with different caps");
        if (other.bins > bins) {
            bins = other.bins;
            counts.resize(bins * STATUSES, 0);
        }
        for (size_t i = 0; i < other.counts.size(); i++)
            counts[i] += other.counts[i];
    }
    void clear() {
        std::fill(counts.begin(), counts.end(), (uint64_t)0);
    }
};

/// §3.5.2's per-shot component structure, over **all** components of the shot's `H`.
///
/// Computed outside every timed window and read by nothing on the decode path. `BallProfile` holds
/// the original and `SpecMatchingProfile` mirrors it, so the campaign accumulator and the
/// row-aligned pybind columns both see one definition.
///
/// There is no resolver classification here any more. §2 deleted every lookup-table resolver: each
/// component is decided by truncated sparse blossom on `H[C]`, so "what the trivial resolver would
/// have done" describes no code that runs. What replaces it is `components_truncated` — the actual
/// per-component status, which is the quantity the escalation predicate is built from.
struct ComponentStats {
    /// 1 iff this shot's components were analysed, so a campaign that left the flag off reports no
    /// component structure rather than a corpus of zeros.
    int measured{0};

    int num_components{0};
    /// Components whose own truncated solve did not finish by `T`. The shot escalates iff this is
    /// non-zero (§2.6), and it is read off the decode's own per-component statuses rather than
    /// recomputed here.
    int components_truncated{0};

    int num_singleton_components{0};
    int num_pair_components{0};
    /// Components of size `>= 3`. Named for the structure rather than for a resolver's reach.
    int num_components_size_ge3{0};
    /// Size `<= 2`, and the defects in them. Structural size classes, kept so the §C series stays
    /// comparable with the campaigns that were run before this branch.
    int num_trivial_components{0};
    int defects_in_trivial_components{0};
    int largest_component_size{0};

    /// `H`'s node count, i.e. the post-preamble defects the components partition. The denominator
    /// of every per-shot fraction of them.
    int component_defects{0};

    /// The solve had nothing to do: `H` has no nodes, so there is no component and no instance to
    /// build. Counted rather than derived from a mean.
    inline bool solver_set_empty() const {
        return measured != 0 && component_defects == 0;
    }
};

/// Per-shot profile. Filled only when a non-null pointer is passed down the decode path, so the
/// hot path pays nothing when profiling is off.
struct SpecMatchingProfile {
    /// Phase 1 on `H` (or on `G`, in the oracle configuration): everything up to and excluding the
    /// harvest/extraction. `BallProfile` splits this further.
    long long phase1_ns{0};
    /// Harvest on the completed path, or the extract-only reduction once §M3.4's bypass is live.
    long long harvest_ns{0};
    /// §M7.7. The one terminal `max_u Y(u)` scan on shots that completed on `H`, on the stock-on-`H`
    /// path only. Counted inside `phase1_ns`'s window and reported separately, because it is the
    /// only cost the certificate adds and the §M7.8 read has to net it off the gating overhead the
    /// swap removes.
    long long dual_scan_ns{0};
    /// The whole escalation call on an escalating shot; 0 otherwise. Kept strictly apart from the
    /// Phase-1 stages so the two are never conflated in a mean (§M3.1).
    long long escalation_ns{0};
    /// The stock exact decode *inside* that call, timed on its own. §M3.2's cost table estimated
    /// `stock` as the mean over all shots; the escalation path runs stock on precisely the shots
    /// that escalate, so this is the exact number and it replaces the estimate.
    long long stock_ns{0};
    /// Measured end to end, deliberately *not* the sum of the parts: a gap between them is
    /// unattributed cost, and noticing it is the point.
    long long total_ns{0};
    /// Stock exact decode on the same shot. Benchmark mode only; 0 otherwise.
    long long exact_reference_ns{0};

    int num_defects{0};
    int residual_size{0};
    int num_trees{0};
    int committed_pairs_frozen{0};
    int committed_pairs_tree{0};
    int committed_boundary{0};
    int largest_tree_size{0};
    int exposed_root_blossoms{0};

    /// §M2.9.6. Harvest's four stages, and the counts and depths that the critical-path model is
    /// built from. Filled only when `Harvester::collect_diagnostics` is set; see `BallProfile` for
    /// what each one means.
    long long harvest_enumerate_ns{0};
    long long harvest_reduce_ns{0};
    long long harvest_base_descent_ns{0};
    long long harvest_shatter_ns{0};
    int blossom_formations{0};
    int max_blossom_nesting_depth{0};
    int matched_blossom_shatters{0};
    /// Profiling builds only: harvest's modelled serial depth, and the solve's, so that the two can
    /// be compared in the units the M2.9 exit checkpoint asks for.
    int harvest_dependent_depth{0};
    int solve_dependent_depth{0};
    int solve_events{0};

    /// §3.5.2, mirrored off `BallProfile` after the shot's timed window has closed. All zero, and
    /// `components.measured == 0`, unless `collect_component_stats` is on and Phase 1 ran on `H`.
    ComponentStats components;

    /// §2. The components this shot's `H` decomposed into, and how many of them truncated. Mirrored
    /// off `BallProfile`, and filled on **every** profiled shot rather than only the analysed ones:
    /// they are the decode's own tally, not part of §3.5.2's post-shot decomposition, and
    /// `any_component_truncated` is read from the second of them.
    int components_total{0};
    int components_truncated{0};

    /// `Sum_S y_S <= exact optimum`. Zero on an escalating shot: the truncated dual is discarded
    /// along with the rest of Phase 1's partial result, so there is nothing to certify against.
    pm::total_weight_int dual_sum_at_truncation{0};
    pm::total_weight_int weight_out{0};
    /// §M7.7. `max_u Y(u)` at completion — the quantity the certificate tests. Zero when `H` could
    /// not complete, and on the truncated-`H` path, where there is no certificate.
    pm::total_weight_int max_dual_at_completion{0};

    /// §M7. 1 iff the shot completed on `H` with `max_u Y(u) <= T_int`, i.e. iff `H`'s answer was
    /// kept as certified globally optimal.
    int certified{0};
    /// §M7. 1 iff `H` admitted no perfect matching — the escalation trigger that announces itself,
    /// as opposed to the completing-but-over-`T` one that does not.
    int h_no_perfect_matching{0};
    /// §M7.7 benchmark mode: would the landed truncated-`H` scheme have escalated this same shot?
    /// Replayed on the identical shot so that `q_current - q_this` is measured rather than asserted.
    /// `..._measured` says the replay actually ran, so that a campaign which left it off reports no
    /// `q_current` at all rather than reporting a spurious zero.
    bool truncated_reference_escalates{false};
    bool truncated_reference_measured{false};

    /// §2.6's **one** escalation predicate: some component of `H` did not resolve within `T`. On
    /// the truncated path that is "trees survived at `T`" in that component; on §M7's it is "the
    /// component's certificate did not hold". Same branch, same downstream contract, and it stays
    /// the quantity invariant 12 is written against.
    ///
    /// Named for what it is rather than for the shot, because after §2 there is no shot-level
    /// truncation: `H` is never solved as one problem, so "the shot truncated" would be a statement
    /// about a solve that does not happen.
    bool any_component_truncated{false};
    /// Whether the shot was re-decoded by stock. Debug-asserted equal to `any_component_truncated`,
    /// which is debug invariant 12 and §2.6's "one predicate, in one place".
    bool escalated{false};
    /// The thread was descheduled during this shot. Its timings no longer *include* the off-CPU
    /// time — `HiResTimer` is thread-scoped — but it resumed on cold caches, so it is still an
    /// outlier. Excluded from percentiles and counted in `contaminated_shot_rate`.
    bool contaminated{false};

    void clear() {
        *this = SpecMatchingProfile();
    }

    /// Copies across everything the harvest already knows, so callers do not restate it.
    void fill_from(const HarvestResult& harvest) {
        residual_size = (int)harvest.residual.size();
        num_trees = harvest.num_trees;
        committed_pairs_frozen = harvest.committed_pairs_frozen;
        committed_pairs_tree = harvest.committed_pairs_tree + harvest.committed_pairs_blossom_cycle;
        committed_boundary = harvest.committed_boundary;
        largest_tree_size = harvest.largest_tree_size;
        exposed_root_blossoms = harvest.exposed_root_blossoms;
        dual_sum_at_truncation = harvest.dual_sum_at_truncation;
        // `any_component_truncated` is deliberately **not** set from the residual. It is the
        // decode's own per-component predicate, and the production path never builds a residual at
        // all — reading one here is how the two would drift apart (§2.6).

        harvest_enumerate_ns = harvest.enumerate_ns;
        harvest_reduce_ns = harvest.reduce_ns;
        harvest_base_descent_ns = harvest.base_descent_ns;
        harvest_shatter_ns = harvest.shatter_ns;
        max_blossom_nesting_depth = harvest.max_blossom_nesting_depth;
        matched_blossom_shatters = harvest.matched_blossom_shatters;
        harvest_dependent_depth = harvest.harvest_dependent_depth;
    }
};

/// Rates need an accumulator that the per-shot struct cannot provide. Accumulated inside
/// `SpecMatchingDecoder::decode_batch` whenever profiling is on (§M6.2).
struct SpecMatchingAggregateStats {
    uint64_t shots{0};
    /// §2.6. **The** escalation count: shots with at least one `TRUNCATED` component. There is no
    /// separate `shots_truncated` any more — it counted a shot-level solve that no longer happens,
    /// and every §M6.2 quantity that used to read it reads this instead.
    uint64_t shots_escalated{0};
    uint64_t shots_zero_defects{0};
    /// §M7. Shots the certificate kept, and the split of the escalating ones by trigger.
    uint64_t shots_certified{0};
    uint64_t shots_h_no_perfect_matching{0};
    /// §M7.7 benchmark mode: shots the **landed truncated scheme** would have escalated, replayed
    /// on the identical corpus. `q_current_on_same_corpus` is this over `shots`, and the design's
    /// `q_this <= q_current` claim is read off the pair rather than assumed.
    uint64_t shots_escalated_truncated_reference{0};
    /// Shots on which that replay actually ran, so a campaign that left `measure_truncated_reference`
    /// off reports no `q_current` rather than reporting zero.
    uint64_t shots_with_truncated_reference{0};
    /// Shots the scheduler interfered with. Never dropped silently: they are counted here and
    /// excluded from the per-shot vectors, and `contaminated_shot_rate` is reported beside every
    /// percentile (§M6 exit checkpoint).
    uint64_t shots_contaminated{0};

    long long sum_phase1_ns{0};
    long long sum_harvest_ns{0};
    /// §M7.7. Inside `sum_phase1_ns`, and reported beside it.
    long long sum_dual_scan_ns{0};
    /// Escalating shots only.
    long long sum_escalation_ns{0};
    long long sum_total_ns{0};
    /// §2.6. End-to-end time of the escalating shots, which is what `sum_total_ns_truncated` used
    /// to be under the other name — the two predicates are one now.
    long long sum_total_ns_escalated{0};
    long long sum_exact_reference_ns{0};
    /// §M3.2: the stock decode measured on precisely the shots that escalate, which replaces that
    /// section's mean-over-all-shots estimate.
    long long sum_stock_ns_on_escalated{0};

    /// Benchmark mode only: kept so percentiles can be computed, and holding uncontaminated shots
    /// only. Long campaigns leave these empty.
    bool keep_per_shot{false};
    std::vector<long long> per_shot_total_ns;
    std::vector<long long> per_shot_exact_ns;
    /// The escalating shots' end-to-end times, kept separately because they are the tail the
    /// streaming budget is sized from and there are too few of them to find by percentile.
    std::vector<long long> per_shot_escalated_total_ns;
    long long max_total_ns{0};
    long long max_escalated_total_ns{0};

    /// Residual sizes, capped, with the last bin acting as the overflow bin.
    static constexpr size_t RESIDUAL_HIST_BINS = 33;
    std::vector<uint64_t> residual_size_hist = std::vector<uint64_t>(RESIDUAL_HIST_BINS, 0);
    uint64_t max_residual_size{0};
    uint64_t max_largest_tree_size{0};

    uint64_t sum_residual_size{0};
    uint64_t sum_num_defects{0};
    uint64_t sum_committed_pairs_frozen{0};
    uint64_t sum_committed_pairs_tree{0};
    uint64_t sum_committed_boundary{0};
    uint64_t sum_exposed_root_blossoms{0};
    uint64_t shots_with_exposed_root_blossoms{0};

    /// §2.6. The components decided across the campaign and how many of them truncated, summed over
    /// **every** shot rather than only the analysed ones: they are the decode's own tally.
    /// `components_truncated / components_total` is the per-component escalation rate, which is a
    /// different and finer quantity than `q`, and `truncated_components_per_escalated_shot` is read
    /// off the pair below (§3.5.1).
    uint64_t components_total{0};
    uint64_t components_truncated{0};
    /// Truncated components summed over the escalating shots alone. Equal to `components_truncated`
    /// by construction — a shot escalates iff it has one — and kept separately so the mean of
    /// §3.5.1 has its own numerator rather than an identity a refactor could quietly break.
    uint64_t sum_truncated_components_on_escalated{0};

    /// §2.6/§3.5.2's joint table, over the components of every analysed shot. Configured once per
    /// campaign, beside the histogram, and defaulted to the same bins the histogram defaults to.
    ComponentStatusTable size_x_status = ComponentStatusTable::with_bins(ComponentHistograms::SIZE_HIST_BINS);

    /// §3.5.2. Running sums over the shots that were analysed, and the distributions.
    ///
    /// `shots_with_component_stats` is the divisor for every mean below, and is *not* `shots`: a
    /// campaign may analyse a subset, and a campaign that left `collect_component_stats` off must
    /// report no component structure rather than a mean over zeros.
    uint64_t shots_with_component_stats{0};
    /// Analysed shots that had at least one `H` node. The divisor of `solver_set_empty_rate`, and
    /// not `shots_with_component_stats`: a shot with no defects has nothing to solve for a reason
    /// that has nothing to do with `H`'s structure, and pooling the two would report the corpus's
    /// zero-defect rate as if it were a property of the decomposition.
    uint64_t shots_with_component_defects{0};
    /// Analysed shots with no `H` node at all, so there was no component and no instance to build.
    uint64_t shots_solver_set_empty{0};
    uint64_t sum_num_components{0};
    uint64_t sum_trivial_components{0};
    uint64_t sum_singleton_components{0};
    uint64_t sum_pair_components{0};
    /// `H`'s node count summed over the analysed shots — the post-preamble defects, which is what
    /// the components actually partition, and not the raw detection-event count `sum_num_defects`.
    uint64_t sum_component_defects{0};
    uint64_t sum_defects_in_trivial_components{0};
    uint64_t sum_components_size_ge3{0};
    /// Per-shot maxima, summed and maximised: the mean says what a typical shot's largest component
    /// looks like, the max says what the largest shot of the campaign held.
    uint64_t sum_largest_component_size{0};
    uint64_t max_component_size{0};
    ComponentHistograms component_hist;

    void reset() {
        *this = SpecMatchingAggregateStats();
    }

    /// Moves the component-size cap onto the histogram and the joint table together, so a campaign
    /// cannot end up with a `size_hist` and a `size_x_status` binned differently.
    void configure_component_tables(size_t size_cap) {
        component_hist.configure(size_cap);
        size_x_status.configure(size_cap + 1);
    }

    /// Removes it from both, the same way and for the same reason. See
    /// `ComponentHistograms::configure_uncapped`.
    void configure_component_tables_uncapped() {
        component_hist.configure_uncapped();
        size_x_status.configure_uncapped();
    }

    /// §3.5.2's distributions for one shot, folded in. Called beside `accumulate` from the profiling
    /// branch of `decode_batch`, and by any driver that decodes shot by shot instead.
    void accumulate_component_histograms(const ComponentHistograms& histograms) {
        component_hist.add(histograms);
    }

    /// §2.6/§3.5.2's joint table for one shot, folded in beside the histograms.
    void accumulate_component_status_table(const ComponentStatusTable& size) {
        size_x_status.add(size);
    }

    void accumulate(const SpecMatchingProfile& profile) {
        shots++;
        // §2.6's one predicate, counted once. `escalated` and `any_component_truncated` are the same
        // fact and are debug-asserted equal where the branch is taken; this counts the branch.
        assert(profile.escalated == profile.any_component_truncated && "§2.6: two escalation predicates disagree");
        if (profile.escalated)
            shots_escalated++;
        if (profile.num_defects == 0)
            shots_zero_defects++;
        if (profile.contaminated)
            shots_contaminated++;
        shots_certified += (uint64_t)profile.certified;
        shots_h_no_perfect_matching += (uint64_t)profile.h_no_perfect_matching;
        if (profile.truncated_reference_measured) {
            shots_with_truncated_reference++;
            if (profile.truncated_reference_escalates)
                shots_escalated_truncated_reference++;
        }

        sum_phase1_ns += profile.phase1_ns;
        sum_harvest_ns += profile.harvest_ns;
        sum_dual_scan_ns += profile.dual_scan_ns;
        sum_total_ns += profile.total_ns;
        sum_exact_reference_ns += profile.exact_reference_ns;
        if (profile.escalated) {
            sum_escalation_ns += profile.escalation_ns;
            sum_stock_ns_on_escalated += profile.stock_ns;
            sum_total_ns_escalated += profile.total_ns;
            sum_truncated_components_on_escalated += (uint64_t)profile.components_truncated;
        }

        if (!profile.contaminated) {
            max_total_ns = std::max(max_total_ns, profile.total_ns);
            if (profile.escalated)
                max_escalated_total_ns = std::max(max_escalated_total_ns, profile.total_ns);
            if (keep_per_shot) {
                per_shot_total_ns.push_back(profile.total_ns);
                per_shot_exact_ns.push_back(profile.exact_reference_ns);
                if (profile.escalated)
                    per_shot_escalated_total_ns.push_back(profile.total_ns);
            }
        }

        size_t bin = std::min((size_t)profile.residual_size, RESIDUAL_HIST_BINS - 1);
        residual_size_hist[bin]++;
        max_residual_size = std::max(max_residual_size, (uint64_t)profile.residual_size);
        max_largest_tree_size = std::max(max_largest_tree_size, (uint64_t)profile.largest_tree_size);

        sum_residual_size += (uint64_t)profile.residual_size;
        sum_num_defects += (uint64_t)profile.num_defects;
        sum_committed_pairs_frozen += (uint64_t)profile.committed_pairs_frozen;
        sum_committed_pairs_tree += (uint64_t)profile.committed_pairs_tree;
        sum_committed_boundary += (uint64_t)profile.committed_boundary;
        sum_exposed_root_blossoms += (uint64_t)profile.exposed_root_blossoms;
        if (profile.exposed_root_blossoms > 0)
            shots_with_exposed_root_blossoms++;

        // §2.6. The decode's own per-component tally, over every shot: the escalation predicate is
        // built from it, so it is not conditioned on the §3.5.2 analysis having run.
        components_total += (uint64_t)profile.components_total;
        components_truncated += (uint64_t)profile.components_truncated;
        assert(
            (profile.components_truncated > 0) == profile.any_component_truncated &&
            "§2.6: the component tally and the escalation predicate disagree");

        const ComponentStats& components = profile.components;
        if (components.measured) {
            shots_with_component_stats++;
            if (components.component_defects > 0)
                shots_with_component_defects++;
            if (components.solver_set_empty())
                shots_solver_set_empty++;
            sum_num_components += (uint64_t)components.num_components;
            sum_trivial_components += (uint64_t)components.num_trivial_components;
            sum_singleton_components += (uint64_t)components.num_singleton_components;
            sum_pair_components += (uint64_t)components.num_pair_components;
            sum_components_size_ge3 += (uint64_t)components.num_components_size_ge3;
            sum_defects_in_trivial_components += (uint64_t)components.defects_in_trivial_components;
            sum_component_defects += (uint64_t)components.component_defects;
            sum_largest_component_size += (uint64_t)components.largest_component_size;
            max_component_size = std::max(max_component_size, (uint64_t)components.largest_component_size);
        }
    }
};

/// The derived quantities of §M6.2's table, computed **once** here rather than re-derived in every
/// benchmark script. `specmatching.summarize()` (§M6.3) wraps this.
struct SpecMatchingSummary {
    /// Escalation rate: how often the shot is re-decoded on `G`. **Always report `shots` beside
    /// it** — a zero below `1 / shots` is a resolution floor, not a measurement (§M3.0). On §M7's
    /// path this is that section's `q_this`, and the bindings expose it under both names.
    double q{0};
    uint64_t shots{0};
    uint64_t shots_escalated{0};
    /// §M7.7. The landed truncated scheme's decision replayed on the identical shots, so that
    /// `q_this <= q_current` is a measurement rather than an assertion. `-1` when the replay was not
    /// run, which is not the same statement as `0`.
    double q_current_on_same_corpus{-1};
    uint64_t shots_with_truncated_reference{0};
    /// §M7. Fraction of shots the certificate kept, and the share of *all* shots that escalated
    /// because `H` had no perfect matching (as against completing with a dual over `T`). The two
    /// escalation triggers are reported apart because only the second one is silent, and it is the
    /// one invariant 4 guards.
    double certified_rate{0};
    double h_no_perfect_matching_rate{0};
    /// §M7.7: the certificate's own cost, per shot. Compare against `c_phase1`.
    double mean_dual_scan_ns{0};
    /// Common-case cost, in ns/shot.
    double c_phase1{0};
    /// Cost conditional on escalating, in ns/shot.
    double c_escalation{0};
    /// `C_phase1 + q * C_escalation`.
    double amortised_mean_ns{0};
    /// `sum_total_ns / shots`. Must reconcile with `amortised_mean_ns`; a gap is unattributed cost.
    double measured_mean_ns{0};
    /// `|amortised - measured| / measured`, so the reconciliation is a number rather than a
    /// judgement call.
    double amortisation_gap{0};
    /// `q * (C_escalation / C_phase1)`: what the fallback costs the mean, as a fraction.
    double amortised_penalty{0};
    /// Expected `≈ 1 + crit_speedup` (§M3.2).
    double escalation_cost_ratio{0};
    /// §M3.2's `stock` term, measured rather than estimated.
    double mean_stock_ns_on_escalated{0};
    /// Headline number, and its tail. At `q ~ 3e-4` the escalation spike sits near p99.97, so p99
    /// and p999 will not show it — hence p9999 and the max.
    double speedup_vs_stock{0};
    double p50_total_ns{0};
    double p99_total_ns{0};
    double p999_total_ns{0};
    double p9999_total_ns{0};
    double max_total_ns{0};
    /// Worst-case escalation latency, which is what a streaming budget is actually sized from.
    double max_escalated_total_ns{0};
    /// Reported beside every percentile above (§M6 exit checkpoint).
    double contaminated_shot_rate{0};

    double mean_residual_size{0};
    double mean_residual_density{0};
    double exposed_root_blossom_rate{0};

    /// §2.6/§3.5.1's per-component escalation reading, over **every** shot: the components decided,
    /// how many truncated, and the mean number of truncated components on an escalating shot. `q`
    /// says how often a shot escalates; these say how much of `H` was responsible.
    uint64_t components_total{0};
    uint64_t components_truncated{0};
    double truncated_component_rate{0};
    double mean_components_per_shot{0};
    double truncated_components_per_escalated_shot{0};

    /// §3.5.2's component structure, derived once here rather than in every benchmark script.
    ///
    /// `shots_with_component_stats` is the divisor of every mean and rate below, and is reported
    /// beside them for the same reason `shots` is reported beside `q`: zero of them means the
    /// campaign did not measure this, which is a different statement from "it measured zero".
    uint64_t shots_with_component_stats{0};
    /// The subset of those that had at least one defect — the divisor of `solver_set_empty_rate`.
    uint64_t shots_with_component_defects{0};
    double mean_components{0};
    double mean_singleton_components{0};
    double mean_pair_components{0};
    double mean_components_size_ge3{0};
    double mean_largest_component_size{0};
    double max_component_size{0};
    /// How much of the defect set sits in a component of size `<= 2`, over `H`'s nodes summed across
    /// the analysed shots. Purely structural: after §2 every defect goes to the solver, so there is
    /// no companion "fraction the solver would still have had to take" — it is 1 by construction.
    double frac_defects_in_trivial_components{0};
    /// Fraction of analysed shots **with at least one defect** on which there was nothing to solve.
    /// Conditioned on having defects so that a corpus's zero-defect rate is not read as a win;
    /// `shots_with_component_defects` is the divisor.
    double solver_set_empty_rate{0};

    /// §3.5.2's distributions, copied out of the accumulator so `summarize` is the one place a
    /// consumer has to look. See `ComponentHistograms` for the binning; `weight_hist_bins_per_T` is
    /// carried alongside so a plot can label the axis without restating the convention.
    uint64_t weight_hist_bins_per_T{ComponentHistograms::BINS_PER_T};
    std::vector<uint64_t> component_size_hist;
    std::vector<uint64_t> h_edge_weight_hist;
    std::vector<uint64_t> boundary_cost_hist;
    /// §2.6's `size x status` table. Row-major `[bin][status]` with `status` 0 = `COMPLETE`,
    /// 1 = `TRUNCATED`, last bin overflowing unless the campaign configured it uncapped.
    ComponentStatusTable size_x_status;
};

inline double percentile_of(std::vector<long long> values, double fraction) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    size_t index = (size_t)(fraction * (double)(values.size() - 1));
    return (double)values[std::min(index, values.size() - 1)];
}

inline SpecMatchingSummary summarize(const SpecMatchingAggregateStats& stats) {
    SpecMatchingSummary summary;
    if (stats.shots == 0)
        return summary;
    double shots = (double)stats.shots;
    summary.shots = stats.shots;
    summary.shots_escalated = stats.shots_escalated;
    summary.q = (double)stats.shots_escalated / shots;
    summary.certified_rate = (double)stats.shots_certified / shots;
    summary.h_no_perfect_matching_rate = (double)stats.shots_h_no_perfect_matching / shots;
    summary.mean_dual_scan_ns = (double)stats.sum_dual_scan_ns / shots;
    summary.shots_with_truncated_reference = stats.shots_with_truncated_reference;
    if (stats.shots_with_truncated_reference) {
        summary.q_current_on_same_corpus =
            (double)stats.shots_escalated_truncated_reference / (double)stats.shots_with_truncated_reference;
    }
    summary.c_phase1 = (double)(stats.sum_phase1_ns + stats.sum_harvest_ns) / shots;
    summary.c_escalation =
        stats.shots_escalated ? (double)stats.sum_escalation_ns / (double)stats.shots_escalated : 0.0;
    summary.amortised_mean_ns = summary.c_phase1 + summary.q * summary.c_escalation;
    summary.measured_mean_ns = (double)stats.sum_total_ns / shots;
    if (summary.measured_mean_ns > 0) {
        summary.amortisation_gap =
            std::abs(summary.amortised_mean_ns - summary.measured_mean_ns) / summary.measured_mean_ns;
    }
    if (summary.c_phase1 > 0)
        summary.amortised_penalty = summary.q * (summary.c_escalation / summary.c_phase1);
    // §2.6. Reads the escalation predicate, not a second truncation counter: the two were the same
    // number under different names, and keeping only one is the point of that section.
    if (stats.shots_escalated && stats.sum_total_ns) {
        summary.escalation_cost_ratio =
            ((double)stats.sum_total_ns_escalated / (double)stats.shots_escalated) / summary.measured_mean_ns;
    }
    if (stats.shots_escalated)
        summary.mean_stock_ns_on_escalated = (double)stats.sum_stock_ns_on_escalated / (double)stats.shots_escalated;
    if (stats.sum_total_ns)
        summary.speedup_vs_stock = (double)stats.sum_exact_reference_ns / (double)stats.sum_total_ns;
    summary.p50_total_ns = percentile_of(stats.per_shot_total_ns, 0.5);
    summary.p99_total_ns = percentile_of(stats.per_shot_total_ns, 0.99);
    summary.p999_total_ns = percentile_of(stats.per_shot_total_ns, 0.999);
    summary.p9999_total_ns = percentile_of(stats.per_shot_total_ns, 0.9999);
    summary.max_total_ns = (double)stats.max_total_ns;
    summary.max_escalated_total_ns = (double)stats.max_escalated_total_ns;
    summary.contaminated_shot_rate = (double)stats.shots_contaminated / shots;
    summary.mean_residual_size = (double)stats.sum_residual_size / shots;
    if (stats.sum_num_defects)
        summary.mean_residual_density = (double)stats.sum_residual_size / (double)stats.sum_num_defects;
    summary.exposed_root_blossom_rate = (double)stats.shots_with_exposed_root_blossoms / shots;

    // §2.6/§3.5.1. The per-component reading, over every shot: the decode tallies its components
    // whether or not the §3.5.2 analysis ran, because the escalation predicate is built from them.
    summary.components_total = stats.components_total;
    summary.components_truncated = stats.components_truncated;
    summary.mean_components_per_shot = (double)stats.components_total / shots;
    if (stats.components_total)
        summary.truncated_component_rate = (double)stats.components_truncated / (double)stats.components_total;
    if (stats.shots_escalated) {
        summary.truncated_components_per_escalated_shot =
            (double)stats.sum_truncated_components_on_escalated / (double)stats.shots_escalated;
    }

    // §3.5.2. Every mean here is over the *analysed* shots, and the histograms are copied out
    // whether or not any shot was analysed, so a consumer always finds the same keys with the same
    // shapes and reads "not measured" off `shots_with_component_stats` rather than off their absence.
    summary.shots_with_component_stats = stats.shots_with_component_stats;
    summary.shots_with_component_defects = stats.shots_with_component_defects;
    if (stats.shots_with_component_defects) {
        summary.solver_set_empty_rate =
            (double)stats.shots_solver_set_empty / (double)stats.shots_with_component_defects;
    }
    summary.max_component_size = (double)stats.max_component_size;
    if (stats.shots_with_component_stats) {
        double analysed = (double)stats.shots_with_component_stats;
        summary.mean_components = (double)stats.sum_num_components / analysed;
        summary.mean_singleton_components = (double)stats.sum_singleton_components / analysed;
        summary.mean_pair_components = (double)stats.sum_pair_components / analysed;
        summary.mean_components_size_ge3 = (double)stats.sum_components_size_ge3 / analysed;
        // A per-shot maximum, so its mean is "the typical shot's largest", not a mean over
        // components — the histogram beside it is what describes the population.
        summary.mean_largest_component_size = (double)stats.sum_largest_component_size / analysed;
    }
    if (stats.sum_component_defects) {
        double defects = (double)stats.sum_component_defects;
        summary.frac_defects_in_trivial_components = (double)stats.sum_defects_in_trivial_components / defects;
    }
    summary.component_size_hist = stats.component_hist.size_hist;
    summary.h_edge_weight_hist = stats.component_hist.edge_weight_hist;
    summary.boundary_cost_hist = stats.component_hist.bcost_hist;
    summary.size_x_status = stats.size_x_status;
    return summary;
}

}  // namespace spec_matching
}  // namespace pm

#endif  // SPECMATCHING_SPEC_MATCHING_PERF_SPEC_MATCHING_PROFILE_H
