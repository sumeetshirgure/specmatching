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
#include <cstdint>
#include <vector>

#include "pyrematching/sparse_blossom/ints.h"
#include "pyrematching/two_phase/truncation/harvest.h"

#if defined(PYREMATCHING_USE_RDTSC)
#include <x86intrin.h>
#endif

namespace pm {
namespace two_phase {

/// Monotonic nanosecond clock for the per-shot profile.
///
/// The default backend is `std::chrono::high_resolution_clock`. Building with
/// `-DPYREMATCHING_USE_RDTSC=ON` switches to a calibrated `rdtsc`, which is cheaper per reading but
/// only meaningful on an invariant-TSC machine; the calibration happens once, lazily.
struct HiResTimer {
#if defined(PYREMATCHING_USE_RDTSC)
    static double ns_per_tick() {
        static const double calibrated = [] {
            auto chrono_start = std::chrono::high_resolution_clock::now();
            uint64_t tsc_start = __rdtsc();
            // Busy-wait rather than sleep: this runs once, and sleeping would measure the
            // scheduler rather than the clock.
            while (std::chrono::high_resolution_clock::now() - chrono_start < std::chrono::milliseconds(20)) {
            }
            uint64_t tsc_end = __rdtsc();
            auto chrono_end = std::chrono::high_resolution_clock::now();
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
    std::chrono::high_resolution_clock::time_point start_time{};

    inline void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    inline long long elapsed_ns() const {
        return (long long)std::chrono::duration_cast<std::chrono::nanoseconds>(
                   std::chrono::high_resolution_clock::now() - start_time)
            .count();
    }
#endif
};

/// Per-shot profile. Filled only when a non-null pointer is passed down the decode path, so the
/// hot path pays nothing when profiling is off.
struct TwoPhaseProfile {
    long long phase1_ns{0};
    long long harvest_ns{0};
    long long inject_ns{0};
    long long inner_blossom_ns{0};
    long long lift_ns{0};
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
    int portal_collisions{0};
    int inner_detection_events{0};

    /// `Sum_S y_S <= exact optimum`.
    pm::total_weight_int dual_sum_at_truncation{0};
    pm::total_weight_int weight_out{0};

    /// Trees survived at `T`, equivalently the residual is non-empty.
    bool truncated{false};
    /// Whether Phase 2 ran. Debug-asserted equal to `truncated`.
    bool phase2_ran{false};
    bool escalated{false};

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
    }
};

/// Rates need an accumulator that the per-shot struct cannot provide. Accumulated over a campaign
/// (inside `decode_batch` once the M3 driver exists; by the M1 exit-artifact harness until then).
struct TwoPhaseAggregateStats {
    uint64_t shots{0};
    uint64_t shots_truncated{0};
    uint64_t shots_escalated{0};
    uint64_t shots_zero_defects{0};

    long long sum_phase1_ns{0};
    long long sum_harvest_ns{0};
    /// inject + inner + lift, over truncated shots only.
    long long sum_phase2_ns{0};
    long long sum_total_ns{0};
    long long sum_total_ns_truncated{0};
    long long sum_exact_reference_ns{0};

    /// Benchmark mode only: kept so percentiles can be computed. Long campaigns leave these empty.
    bool keep_per_shot{false};
    std::vector<long long> per_shot_total_ns;
    std::vector<long long> per_shot_exact_ns;

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

        sum_phase1_ns += profile.phase1_ns;
        sum_harvest_ns += profile.harvest_ns;
        sum_total_ns += profile.total_ns;
        sum_exact_reference_ns += profile.exact_reference_ns;
        if (profile.truncated) {
            sum_phase2_ns += profile.inject_ns + profile.inner_blossom_ns + profile.lift_ns;
            sum_total_ns_truncated += profile.total_ns;
        }

        if (keep_per_shot) {
            per_shot_total_ns.push_back(profile.total_ns);
            per_shot_exact_ns.push_back(profile.exact_reference_ns);
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

/// The derived quantities of the design's summary table, computed **once** here rather than
/// re-derived in every benchmark script. `pyrematching.summarize()` (M5.3) wraps this.
struct TwoPhaseSummary {
    /// Fallback rate: how often Phase 2 has to run.
    double q{0};
    /// Common-case cost, in ns/shot. Expected to be ~ the stock blossom cost.
    double c_phase1{0};
    /// Cost conditional on falling back, in ns/shot.
    double c_phase2{0};
    /// `C_phase1 + q * C_phase2`.
    double amortised_mean_ns{0};
    /// `sum_total_ns / shots`. Must reconcile with `amortised_mean_ns`; a gap is unattributed cost.
    double measured_mean_ns{0};
    /// The tail heaviness this design targets.
    double fallback_cost_ratio{0};
    /// Headline number, and its tail.
    double speedup_vs_stock{0};
    double p50_total_ns{0};
    double p99_total_ns{0};
    double p999_total_ns{0};
    double p_escalated{0};
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
    summary.q = (double)stats.shots_truncated / shots;
    summary.c_phase1 = (double)(stats.sum_phase1_ns + stats.sum_harvest_ns) / shots;
    summary.c_phase2 = stats.shots_truncated ? (double)stats.sum_phase2_ns / (double)stats.shots_truncated : 0.0;
    summary.amortised_mean_ns = summary.c_phase1 + summary.q * summary.c_phase2;
    summary.measured_mean_ns = (double)stats.sum_total_ns / shots;
    if (stats.shots_truncated && stats.sum_total_ns) {
        summary.fallback_cost_ratio =
            ((double)stats.sum_total_ns_truncated / (double)stats.shots_truncated) / summary.measured_mean_ns;
    }
    if (stats.sum_total_ns)
        summary.speedup_vs_stock = (double)stats.sum_exact_reference_ns / (double)stats.sum_total_ns;
    summary.p50_total_ns = percentile_of(stats.per_shot_total_ns, 0.5);
    summary.p99_total_ns = percentile_of(stats.per_shot_total_ns, 0.99);
    summary.p999_total_ns = percentile_of(stats.per_shot_total_ns, 0.999);
    summary.p_escalated = (double)stats.shots_escalated / shots;
    summary.mean_residual_size = (double)stats.sum_residual_size / shots;
    if (stats.sum_num_defects)
        summary.mean_residual_density = (double)stats.sum_residual_size / (double)stats.sum_num_defects;
    summary.exposed_root_blossom_rate = (double)stats.shots_with_exposed_root_blossoms / shots;
    return summary;
}

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_PERF_TWO_PHASE_PROFILE_H
