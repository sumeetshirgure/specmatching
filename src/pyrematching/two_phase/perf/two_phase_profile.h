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

#ifndef PYREMATCHING_TWO_PHASE_PERF_TWO_PHASE_PROFILE_H
#define PYREMATCHING_TWO_PHASE_PERF_TWO_PHASE_PROFILE_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "pyrematching/perf/thread_timer.h"
#include "pyrematching/sparse_blossom/ints.h"
#include "pyrematching/two_phase/truncation/harvest.h"

#if defined(__linux__)
#include <sys/resource.h>
#endif

/// Linux hands out a per-thread context-switch counter and nothing else does, so `RUSAGE_THREAD` is
/// the preferred mechanism where it exists and `PYREMATCHING_PREEMPTION_CLOCK_GAP` is the fallback
/// used everywhere else that has a thread-scoped clock — macOS in particular. See `PreemptionProbe`.
///
/// `-DPYREMATCHING_FORCE_CLOCK_GAP_PROBE` selects the fallback on a machine that has the counters,
/// which is how the macOS path is exercised from a Linux CI box — the mechanism is chosen at compile
/// time, so without it there is no way to run that code at all.
#if defined(PYREMATCHING_FORCE_CLOCK_GAP_PROBE) && defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
#define PYREMATCHING_PREEMPTION_CLOCK_GAP 1
#elif defined(__linux__) && defined(RUSAGE_THREAD)
#define PYREMATCHING_PREEMPTION_SWITCH_COUNTERS 1
#elif defined(PYREMATCHING_HAVE_THREAD_CPUTIME)
#define PYREMATCHING_PREEMPTION_CLOCK_GAP 1
#endif

namespace pm {
namespace two_phase {

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
/// The slack is `PYREMATCHING_PREEMPTION_SLACK_NS`, default 2000. Well above the noise of four clock
/// readings, well below the tens of microseconds a real deschedule costs.
///
/// Run against the switch counters on the same intervals on a Linux box (which is what
/// `PYREMATCHING_FORCE_CLOCK_GAP_PROBE` is for), the gap flags a **superset**: over 200 intervals of
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
#if defined(PYREMATCHING_PREEMPTION_SWITCH_COUNTERS)
        struct rusage usage;
        if (getrusage(RUSAGE_THREAD, &usage) == 0) {
            voluntary = usage.ru_nvcsw;
            involuntary = usage.ru_nivcsw;
            supported = true;
            return;
        }
#elif defined(PYREMATCHING_PREEMPTION_CLOCK_GAP)
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
#if defined(PYREMATCHING_PREEMPTION_SWITCH_COUNTERS)
        return voluntary != before.voluntary || involuntary != before.involuntary;
#elif defined(PYREMATCHING_PREEMPTION_CLOCK_GAP)
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
#if defined(PYREMATCHING_PREEMPTION_SWITCH_COUNTERS)
        return "rusage_thread_switches";
#elif defined(PYREMATCHING_PREEMPTION_CLOCK_GAP)
        return "wall_minus_thread_cpu";
#else
        return "none";
#endif
    }

    /// Nanoseconds of off-CPU time an interval is allowed before the shot is called contaminated.
    /// Meaningless under the switch-counter mechanism, which does not measure a duration.
    static inline uint64_t slack_ns() {
        static const uint64_t slack = [] {
            const char* value = std::getenv("PYREMATCHING_PREEMPTION_SLACK_NS");
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

/// §C.2's distributions, over one shot's components, `H` edges and defects.
///
/// Kept apart from the per-shot profile because their unit is the component, the edge or the node
/// rather than the shot: the profile stays a flat row of scalars the pybind layer can column-ise,
/// and these are folded straight into the campaign accumulator instead. One instance is reused
/// across shots, so a profiled shot allocates nothing for them either.
///
/// **Unit.** The three weight histograms are binned in *sixteenths of the horizon* — bin `k` counts
/// values in `[k * T / 16, (k + 1) * T / 16)`, and the last bin overflows — so all three share one
/// axis and a plot can label it without knowing the normalising constant. `H`'s defect-defect edges
/// are `<= 2T` by construction, which is exactly the last bin of `edge_weight_hist`.
///
/// Filled outside every timed region, and never counted toward any reported latency.
struct ComponentHistograms {
    static constexpr size_t BINS_PER_T = 16;
    /// Component sizes `0..31`, last bin overflowing. Bin 0 is always empty: a component has at
    /// least one member.
    static constexpr size_t SIZE_HIST_BINS = 33;
    /// `H`-node degrees, in the defect-defect adjacency alone.
    static constexpr size_t DEGREE_HIST_BINS = 33;
    /// `w_int` over `H`'s defect-defect edges: `0 .. 2T`, which is the whole range.
    static constexpr size_t WEIGHT_HIST_BINS = 33;
    /// `bcost_int` over defects that have a boundary entry within `R`. Values past `2T` — legal,
    /// since `R >= 2 * T_max` — land in the overflow bin.
    static constexpr size_t BCOST_HIST_BINS = 33;
    /// Weighted `H`-subgraph diameters: `0 .. 4T`, then overflow.
    static constexpr size_t DIAMETER_HIST_BINS = 65;

    std::vector<uint64_t> size_hist = std::vector<uint64_t>(SIZE_HIST_BINS, 0);
    std::vector<uint64_t> degree_hist = std::vector<uint64_t>(DEGREE_HIST_BINS, 0);
    std::vector<uint64_t> edge_weight_hist = std::vector<uint64_t>(WEIGHT_HIST_BINS, 0);
    std::vector<uint64_t> bcost_hist = std::vector<uint64_t>(BCOST_HIST_BINS, 0);
    std::vector<uint64_t> diameter_hist = std::vector<uint64_t>(DIAMETER_HIST_BINS, 0);

    /// Zeroes every bin without touching the buffers, so the next shot reuses the same storage.
    void clear() {
        std::fill(size_hist.begin(), size_hist.end(), (uint64_t)0);
        std::fill(degree_hist.begin(), degree_hist.end(), (uint64_t)0);
        std::fill(edge_weight_hist.begin(), edge_weight_hist.end(), (uint64_t)0);
        std::fill(bcost_hist.begin(), bcost_hist.end(), (uint64_t)0);
        std::fill(diameter_hist.begin(), diameter_hist.end(), (uint64_t)0);
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
        for (size_t i = 0; i < size_hist.size(); i++)
            size_hist[i] += other.size_hist[i];
        for (size_t i = 0; i < degree_hist.size(); i++)
            degree_hist[i] += other.degree_hist[i];
        for (size_t i = 0; i < edge_weight_hist.size(); i++)
            edge_weight_hist[i] += other.edge_weight_hist[i];
        for (size_t i = 0; i < bcost_hist.size(); i++)
            bcost_hist[i] += other.bcost_hist[i];
        for (size_t i = 0; i < diameter_hist.size(); i++)
            diameter_hist[i] += other.diameter_hist[i];
    }
};

/// §C.1's per-shot component structure, over **all** components of the shot's `H`.
///
/// Computed outside every timed window and read by nothing on the decode path. `BallProfile` holds
/// the original and `TwoPhaseProfile` mirrors it, so the campaign accumulator and the row-aligned
/// pybind columns both see one definition.
///
/// The last three fields are the classification of §A.3 — what the trivial resolver *would* commit,
/// leave as residual, and hand to the solver — evaluated for its counts alone. Nothing is committed
/// and no decode output depends on them.
struct ComponentStats {
    /// 1 iff this shot's components were analysed, so a campaign that left the flag off reports no
    /// component structure rather than a corpus of zeros.
    int measured{0};

    int num_components{0};
    int num_trivial_components{0};
    int num_singleton_components{0};
    int num_pair_components{0};
    int num_nontrivial_components{0};
    int defects_in_trivial_components{0};
    int largest_component_size{0};
    /// Weighted `H`-subgraph diameter, in integer time units. See `component_diameter` for what
    /// "confined to the component" means and why it is not a `G` distance.
    int max_component_diameter_wint{0};
    /// Components too large for the diameter to be computed (§C.2's cap).
    int diameter_uncomputed_components{0};
    /// Components with at least one member whose `bcost_int <= T_int`.
    int num_boundary_touching_components{0};

    int defects_committed_trivially{0};
    int defects_residual_trivially{0};
    int defects_to_solver{0};

    /// Every component resolved trivially, i.e. the solver's input would have been empty and both
    /// the `Mwpm(H)` build and the harvest would have been skipped outright (§A.4). This is the
    /// shot class the `blossom_on_h_ns` / `harvest_ns` reduction comes from, so it is counted
    /// rather than derived from a mean.
    inline bool solver_set_empty() const {
        return measured != 0 && defects_to_solver == 0;
    }
};

/// Per-shot profile. Filled only when a non-null pointer is passed down the decode path, so the
/// hot path pays nothing when profiling is off.
struct TwoPhaseProfile {
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

    /// §C.1, mirrored off `BallProfile` after the shot's timed window has closed. All zero, and
    /// `components.measured == 0`, unless `collect_component_stats` is on and Phase 1 ran on `H`.
    ComponentStats components;

    /// §A's run label, mirrored off `BallProfile`: the `k` this shot decoded at, and the defects the
    /// small-component resolver settled off the solver. Labels on the latency fields above, never
    /// latency themselves — the resolve is a serial pre-pass excluded from them by measurement
    /// scope. See `BallProfile` for what `defects_resolved_small` counts.
    int prune_component_max_size{0};
    int defects_resolved_small{0};

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

    /// Phase 1 did not yield a usable answer, so the shot escalates. On the truncated-`H` path that
    /// is "trees survived at `T`", equivalently "the residual is non-empty"; on §M7's path it is
    /// "the certificate did not hold". Same branch, same downstream contract, and it stays the
    /// quantity invariant 12 is written against.
    bool truncated{false};
    /// Whether the shot was re-decoded by stock. Debug-asserted equal to `truncated`, which is
    /// debug invariant 12.
    bool escalated{false};
    /// The thread was descheduled during this shot. Its timings no longer *include* the off-CPU
    /// time — `HiResTimer` is thread-scoped — but it resumed on cold caches, so it is still an
    /// outlier. Excluded from percentiles and counted in `contaminated_shot_rate`.
    bool contaminated{false};

    void clear() {
        *this = TwoPhaseProfile();
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
        truncated = !harvest.residual.empty();

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
/// `TwoPhaseDecoder::decode_batch` whenever profiling is on (§M6.2).
struct TwoPhaseAggregateStats {
    uint64_t shots{0};
    uint64_t shots_truncated{0};
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
    long long sum_total_ns_truncated{0};
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

    /// §C.3. Running sums of §C.1 over the shots that were analysed, and the §C.2 distributions.
    ///
    /// `shots_with_component_stats` is the divisor for every mean below, and is *not* `shots`: a
    /// campaign may analyse a subset, and a campaign that left `collect_component_stats` off must
    /// report no component structure rather than a mean over zeros.
    uint64_t shots_with_component_stats{0};
    /// Analysed shots that had at least one `H` node. The divisor of `solver_set_empty_rate`, and
    /// not `shots_with_component_stats`: a shot with no defects has an empty solver set for a reason
    /// that has nothing to do with the prune, and pooling the two would quote a win that is really
    /// the corpus's zero-defect rate.
    uint64_t shots_with_component_defects{0};
    /// Shots whose components were all resolved trivially, so the solver and the harvest would have
    /// been skipped entirely (§A.4). Zero-defect shots are in here too — the solve really is skipped
    /// on them — which is why the rate above is conditioned rather than taken over every shot.
    uint64_t shots_solver_set_empty{0};
    /// Shots with at least one trivially-residual defect — a singleton whose boundary sits past the
    /// horizon, which forces escalation on its own (§A.5).
    uint64_t shots_with_trivial_residual{0};
    uint64_t sum_num_components{0};
    uint64_t sum_trivial_components{0};
    uint64_t sum_singleton_components{0};
    uint64_t sum_pair_components{0};
    uint64_t sum_nontrivial_components{0};
    uint64_t sum_boundary_touching_components{0};
    uint64_t sum_diameter_uncomputed_components{0};
    /// The denominator of the two fractions §D asks for. `H`'s node count summed over the analysed
    /// shots — the post-preamble defects, which is what the components actually partition — and not
    /// the raw detection-event count `sum_num_defects`.
    uint64_t sum_component_defects{0};
    uint64_t sum_defects_in_trivial_components{0};
    uint64_t sum_defects_committed_trivially{0};
    uint64_t sum_defects_residual_trivially{0};
    uint64_t sum_defects_to_solver{0};
    /// Per-shot maxima, summed and maximised: the mean says what a typical shot's worst component
    /// looks like, the max says what the worst shot of the campaign held.
    uint64_t sum_largest_component_size{0};
    uint64_t max_component_size{0};
    uint64_t sum_max_component_diameter_wint{0};
    uint64_t max_component_diameter_wint{0};
    ComponentHistograms component_hist;

    /// §A's run label. `k` is a property of the campaign rather than of a shot, so it is recorded
    /// rather than summed — the last shot folded in wins, and pooling two campaigns decoded at
    /// different `k` is a category error the reader has to avoid. `-1` until a shot is accumulated,
    /// which is not the same statement as `0` (a run that resolved nothing off the solver).
    int prune_component_max_size{-1};
    /// Defects the resolver settled off the solver, summed over every accumulated shot. Not a
    /// latency: the resolve is a serial pre-pass outside the reported stages by scope.
    uint64_t sum_defects_resolved_small{0};

    void reset() {
        *this = TwoPhaseAggregateStats();
    }

    /// §C.2's distributions for one shot, folded in. Called beside `accumulate` from the profiling
    /// branch of `decode_batch`, and by any driver that decodes shot by shot instead.
    void accumulate_component_histograms(const ComponentHistograms& histograms) {
        component_hist.add(histograms);
    }

    void accumulate(const TwoPhaseProfile& profile) {
        shots++;
        if (profile.truncated)
            shots_truncated++;
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
        }
        if (profile.truncated)
            sum_total_ns_truncated += profile.total_ns;

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

        // §A. `k` labels the campaign; the defect count is the resolver's load over it.
        prune_component_max_size = profile.prune_component_max_size;
        sum_defects_resolved_small += (uint64_t)profile.defects_resolved_small;

        const ComponentStats& components = profile.components;
        if (components.measured) {
            shots_with_component_stats++;
            if (components.defects_committed_trivially + components.defects_residual_trivially +
                    components.defects_to_solver >
                0)
                shots_with_component_defects++;
            if (components.solver_set_empty())
                shots_solver_set_empty++;
            if (components.defects_residual_trivially > 0)
                shots_with_trivial_residual++;
            sum_num_components += (uint64_t)components.num_components;
            sum_trivial_components += (uint64_t)components.num_trivial_components;
            sum_singleton_components += (uint64_t)components.num_singleton_components;
            sum_pair_components += (uint64_t)components.num_pair_components;
            sum_nontrivial_components += (uint64_t)components.num_nontrivial_components;
            sum_boundary_touching_components += (uint64_t)components.num_boundary_touching_components;
            sum_diameter_uncomputed_components += (uint64_t)components.diameter_uncomputed_components;
            sum_defects_in_trivial_components += (uint64_t)components.defects_in_trivial_components;
            sum_defects_committed_trivially += (uint64_t)components.defects_committed_trivially;
            sum_defects_residual_trivially += (uint64_t)components.defects_residual_trivially;
            sum_defects_to_solver += (uint64_t)components.defects_to_solver;
            // The components partition `H`'s nodes, so their three classifications sum to the node
            // count — which is the denominator §D's fractions are taken over.
            sum_component_defects += (uint64_t)components.defects_committed_trivially +
                                     (uint64_t)components.defects_residual_trivially +
                                     (uint64_t)components.defects_to_solver;
            sum_largest_component_size += (uint64_t)components.largest_component_size;
            max_component_size = std::max(max_component_size, (uint64_t)components.largest_component_size);
            sum_max_component_diameter_wint += (uint64_t)components.max_component_diameter_wint;
            max_component_diameter_wint =
                std::max(max_component_diameter_wint, (uint64_t)components.max_component_diameter_wint);
        }
    }
};

/// The derived quantities of §M6.2's table, computed **once** here rather than re-derived in every
/// benchmark script. `pyrematching.summarize()` (§M6.3) wraps this.
struct TwoPhaseSummary {
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

    /// §C/§D's component structure, derived once here rather than in every benchmark script.
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
    double mean_nontrivial_components{0};
    double mean_largest_component_size{0};
    double max_component_size{0};
    double mean_max_component_diameter_wint{0};
    double max_component_diameter_wint{0};
    double boundary_touching_component_fraction{0};
    double diameter_uncomputed_component_fraction{0};
    /// The two fractions §D asks for, over `H`'s nodes summed across the analysed shots: how much of
    /// the defect set sits in a trivial (size `<= 2`) component, and how much of it the solver would
    /// still have had to take. They do not sum to 1 — an ambiguous trivial component is in the first
    /// and in the second, which is exactly the gap between "trivial" and "trivially resolvable".
    double frac_defects_in_trivial_components{0};
    double frac_defects_to_solver{0};
    double frac_defects_committed_trivially{0};
    double frac_defects_residual_trivially{0};
    /// Fraction of analysed shots **with at least one defect** on which the solver's input would
    /// have been empty, so the `Mwpm(H)` build, the solve and the harvest would not have run at all
    /// (§A.4). Conditioned on having defects so that a corpus's zero-defect rate is not read as the
    /// prune's win; `shots_with_component_defects` is the divisor.
    double solver_set_empty_rate{0};
    double trivial_residual_rate{0};

    /// §A's run label. `k` is what the campaign decoded at — the solver saw only components of size
    /// `> k` — and `mean_defects_resolved_small` is how many of `H`'s defects the small-component
    /// resolver settled off the solver per shot, over every shot, sizes `1..k` included.
    ///
    /// Neither is a latency. The resolve is a serial pre-pass on the critical path, excluded from
    /// the latency this summary reports by measurement scope, to be measured separately. `k = -1`
    /// means no shot was accumulated, which is not the same statement as `k = 0`.
    int prune_component_max_size{-1};
    double mean_defects_resolved_small{0};

    /// §C.2's distributions, copied out of the accumulator so `summarize` is the one place a
    /// consumer has to look. See `ComponentHistograms` for the binning; `weight_hist_bins_per_T` is
    /// carried alongside so a plot can label the axis without restating the convention.
    uint64_t weight_hist_bins_per_T{ComponentHistograms::BINS_PER_T};
    std::vector<uint64_t> component_size_hist;
    std::vector<uint64_t> component_diameter_hist;
    std::vector<uint64_t> h_degree_hist;
    std::vector<uint64_t> h_edge_weight_hist;
    std::vector<uint64_t> boundary_cost_hist;
};

inline double percentile_of(std::vector<long long> values, double fraction) {
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    size_t index = (size_t)(fraction * (double)(values.size() - 1));
    return (double)values[std::min(index, values.size() - 1)];
}

inline TwoPhaseSummary summarize(const TwoPhaseAggregateStats& stats) {
    TwoPhaseSummary summary;
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
    if (stats.shots_truncated && stats.sum_total_ns) {
        summary.escalation_cost_ratio =
            ((double)stats.sum_total_ns_truncated / (double)stats.shots_truncated) / summary.measured_mean_ns;
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

    // §C/§D. Every mean here is over the *analysed* shots, and the histograms are copied out
    // whether or not any shot was analysed, so a consumer always finds the same keys with the same
    // shapes and reads "not measured" off `shots_with_component_stats` rather than off their absence.
    summary.shots_with_component_stats = stats.shots_with_component_stats;
    summary.shots_with_component_defects = stats.shots_with_component_defects;
    if (stats.shots_with_component_defects) {
        summary.solver_set_empty_rate =
            (double)stats.shots_solver_set_empty / (double)stats.shots_with_component_defects;
    }
    summary.max_component_size = (double)stats.max_component_size;
    summary.max_component_diameter_wint = (double)stats.max_component_diameter_wint;
    if (stats.shots_with_component_stats) {
        double analysed = (double)stats.shots_with_component_stats;
        summary.mean_components = (double)stats.sum_num_components / analysed;
        summary.mean_singleton_components = (double)stats.sum_singleton_components / analysed;
        summary.mean_pair_components = (double)stats.sum_pair_components / analysed;
        summary.mean_nontrivial_components = (double)stats.sum_nontrivial_components / analysed;
        summary.trivial_residual_rate = (double)stats.shots_with_trivial_residual / analysed;
        // Both are per-shot maxima, so their means are "the typical shot's worst", not a mean over
        // components — the histograms beside them are what describes the population.
        summary.mean_largest_component_size = (double)stats.sum_largest_component_size / analysed;
        summary.mean_max_component_diameter_wint = (double)stats.sum_max_component_diameter_wint / analysed;
    }
    if (stats.sum_num_components) {
        double components = (double)stats.sum_num_components;
        summary.boundary_touching_component_fraction = (double)stats.sum_boundary_touching_components / components;
        summary.diameter_uncomputed_component_fraction = (double)stats.sum_diameter_uncomputed_components / components;
    }
    if (stats.sum_component_defects) {
        double defects = (double)stats.sum_component_defects;
        summary.frac_defects_in_trivial_components = (double)stats.sum_defects_in_trivial_components / defects;
        summary.frac_defects_to_solver = (double)stats.sum_defects_to_solver / defects;
        summary.frac_defects_committed_trivially = (double)stats.sum_defects_committed_trivially / defects;
        summary.frac_defects_residual_trivially = (double)stats.sum_defects_residual_trivially / defects;
    }
    // §A's run label. Over every shot, not only the analysed ones: the resolver ran on all of them,
    // and `collect_component_stats` has nothing to do with whether it did.
    summary.prune_component_max_size = stats.prune_component_max_size;
    summary.mean_defects_resolved_small = (double)stats.sum_defects_resolved_small / shots;
    summary.component_size_hist = stats.component_hist.size_hist;
    summary.component_diameter_hist = stats.component_hist.diameter_hist;
    summary.h_degree_hist = stats.component_hist.degree_hist;
    summary.h_edge_weight_hist = stats.component_hist.edge_weight_hist;
    summary.boundary_cost_hist = stats.component_hist.bcost_hist;
    return summary;
}

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_PERF_TWO_PHASE_PROFILE_H
