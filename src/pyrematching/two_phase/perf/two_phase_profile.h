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
#include <chrono>
#include <cmath>
#include <cstdint>
#include <vector>

#include "pyrematching/sparse_blossom/ints.h"
#include "pyrematching/two_phase/truncation/harvest.h"

#if defined(PYREMATCHING_USE_RDTSC)
#include <x86intrin.h>
#endif

#if defined(__linux__)
#include <sys/resource.h>
#endif

namespace pm {
namespace two_phase {

/// Monotonic nanosecond clock for the per-shot profile.
///
/// The default backend is `std::chrono::steady_clock`, named explicitly rather than through
/// `high_resolution_clock` — the latter is a *typedef* for `system_clock` on libstdc++, which is
/// not monotonic and can step backwards under NTP, producing negative or absurd per-shot times in
/// exactly the tail percentiles this design cares about (§M6.4 timer discipline).
///
/// Building with `-DPYREMATCHING_USE_RDTSC=ON` switches to a calibrated `rdtsc`, which is cheaper
/// per reading but only meaningful on an invariant-TSC machine; the calibration happens once,
/// lazily.
struct HiResTimer {
#if defined(PYREMATCHING_USE_RDTSC)
    static double ns_per_tick() {
        static const double calibrated = [] {
            auto chrono_start = std::chrono::steady_clock::now();
            uint64_t tsc_start = __rdtsc();
            // Busy-wait rather than sleep: this runs once, and sleeping would measure the
            // scheduler rather than the clock.
            while (std::chrono::steady_clock::now() - chrono_start < std::chrono::milliseconds(20)) {
            }
            uint64_t tsc_end = __rdtsc();
            auto chrono_end = std::chrono::steady_clock::now();
            auto elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(chrono_end - chrono_start).count();
            return (double)elapsed_ns / (double)(tsc_end - tsc_start);
        }();
        return calibrated;
    }

    uint64_t start_tick{0};

    inline void start() {
        start_tick = __rdtsc();
    }
    inline long long elapsed_ns() const {
        return (long long)((double)(__rdtsc() - start_tick) * ns_per_tick());
    }
#else
    std::chrono::steady_clock::time_point start_time{};

    inline void start() {
        start_time = std::chrono::steady_clock::now();
    }
    inline long long elapsed_ns() const {
        return (long long)std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::steady_clock::now() - start_time)
            .count();
    }
#endif
};

/// Was this thread descheduled while the shot was being timed? (§M6.4 timer discipline.)
///
/// One preemption defines p999 outright at 2000 shots, so a percentile quoted without this is a
/// measurement of the scheduler. `RUSAGE_THREAD` is Linux-specific; elsewhere the probe reports
/// "not contaminated" and `contaminated_shot_rate` reads 0, which the artifact says explicitly
/// rather than implying the machine was quiet.
struct PreemptionProbe {
    long voluntary{0};
    long involuntary{0};
    bool supported{false};

    inline void sample() {
#if defined(__linux__) && defined(RUSAGE_THREAD)
        struct rusage usage;
        if (getrusage(RUSAGE_THREAD, &usage) == 0) {
            voluntary = usage.ru_nvcsw;
            involuntary = usage.ru_nivcsw;
            supported = true;
            return;
        }
#endif
        supported = false;
    }

    /// True when a context switch happened between `before` and this sample.
    inline bool switched_since(const PreemptionProbe& before) const {
        if (!supported || !before.supported)
            return false;
        return voluntary != before.voluntary || involuntary != before.involuntary;
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
    /// The thread was descheduled during this shot, so its time is the scheduler's, not the
    /// decoder's. Excluded from percentiles and counted in `contaminated_shot_rate`.
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

    void reset() {
        *this = TwoPhaseAggregateStats();
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
    return summary;
}

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_PERF_TWO_PHASE_PROFILE_H
