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

#ifndef PYREMATCHING_TWO_PHASE_PERF_BALL_PROFILE_H
#define PYREMATCHING_TWO_PHASE_PERF_BALL_PROFILE_H

#include <algorithm>
#include <cstdint>
#include <vector>

#include "pyrematching/two_phase/perf/two_phase_profile.h"

namespace pm {
namespace two_phase {

/// Per-shot profile of the ball-graph front end (§M2.7). Filled only when a non-null pointer is
/// passed down, so the hot path pays nothing when profiling is off.
struct BallProfile {
    /// Walking the balls of the shot's defects and testing the syndrome.
    long long intersect_ns{0};
    /// Turning the hits into `H`: canonical edge ordering and the arena bookkeeping.
    long long h_build_ns{0};
    /// Rewriting the `pm::Mwpm` on `H`.
    long long mwpm_build_ns{0};
    /// `process_timeline_until_horizon` on `H`.
    long long blossom_on_h_ns{0};
    long long harvest_ns{0};
    /// Measured end to end, deliberately *not* the sum of the parts: the gap is unattributed cost.
    long long total_ns{0};
    /// M1's Phase 1 + harvest on `G` for the same shot. Benchmark / verification mode only.
    long long g_reference_ns{0};

    int n_defects{0};
    int h_nodes{0};
    int h_edges{0};
    int h_boundary_edges{0};
    int max_degree{0};
    double mean_degree{0};

    /// Shells actually materialised, and rungs of the restart ladder taken (§M2.8). With the
    /// default `shell_width = 0` these are 1 and 0 on every shot.
    int shells_materialized{1};
    int restarts{0};

    /// Shots where the ball graph's observable bytes differed from `G`'s. Level 2 of §M2.6 counts
    /// these and certifies each one homologically trivial.
    int mask_divergences{0};

    /// Ties resolved differently from `G`: the same optimum reached by a different route. Filled
    /// only when `verify_against_g` is on. See `BallDecoder::verify_level1` for why these are
    /// recorded rather than failed, and what *is* still a hard failure.
    int residual_ties{0};
    int pairing_ties{0};

    void clear() {
        *this = BallProfile();
    }
};

/// Campaign-level accumulator, mirroring `TwoPhaseAggregateStats`. `speedup_vs_m1` is the headline
/// number of the M2 exit read.
struct BallAggregateStats {
    uint64_t shots{0};
    uint64_t shots_truncated{0};
    uint64_t shots_zero_defects{0};

    long long sum_intersect_ns{0};
    long long sum_h_build_ns{0};
    long long sum_mwpm_build_ns{0};
    long long sum_blossom_on_h_ns{0};
    long long sum_harvest_ns{0};
    long long sum_total_ns{0};
    long long sum_g_reference_ns{0};

    uint64_t sum_n_defects{0};
    uint64_t sum_h_nodes{0};
    uint64_t sum_h_edges{0};
    uint64_t sum_h_boundary_edges{0};
    uint64_t max_degree{0};
    uint64_t sum_residual_size{0};
    uint64_t sum_restarts{0};
    uint64_t mask_divergences{0};
    uint64_t residual_ties{0};
    uint64_t pairing_ties{0};

    bool keep_per_shot{false};
    std::vector<long long> per_shot_total_ns;
    std::vector<long long> per_shot_g_reference_ns;

    void reset() {
        *this = BallAggregateStats();
    }

    void accumulate(const BallProfile& profile, const HarvestResult& harvest) {
        shots++;
        if (!harvest.residual.empty())
            shots_truncated++;
        if (profile.n_defects == 0)
            shots_zero_defects++;

        sum_intersect_ns += profile.intersect_ns;
        sum_h_build_ns += profile.h_build_ns;
        sum_mwpm_build_ns += profile.mwpm_build_ns;
        sum_blossom_on_h_ns += profile.blossom_on_h_ns;
        sum_harvest_ns += profile.harvest_ns;
        sum_total_ns += profile.total_ns;
        sum_g_reference_ns += profile.g_reference_ns;

        sum_n_defects += (uint64_t)profile.n_defects;
        sum_h_nodes += (uint64_t)profile.h_nodes;
        sum_h_edges += (uint64_t)profile.h_edges;
        sum_h_boundary_edges += (uint64_t)profile.h_boundary_edges;
        max_degree = std::max(max_degree, (uint64_t)profile.max_degree);
        sum_residual_size += harvest.residual.size();
        sum_restarts += (uint64_t)profile.restarts;
        mask_divergences += (uint64_t)profile.mask_divergences;
        residual_ties += (uint64_t)profile.residual_ties;
        pairing_ties += (uint64_t)profile.pairing_ties;

        if (keep_per_shot) {
            per_shot_total_ns.push_back(profile.total_ns);
            per_shot_g_reference_ns.push_back(profile.g_reference_ns);
        }
    }
};

/// The derived quantities of §M2.7, computed once here rather than in every benchmark script.
struct BallSummary {
    double speedup_vs_m1{0};
    double p50_total_ns{0};
    double p99_total_ns{0};
    double p999_total_ns{0};
    double mean_total_ns{0};
    double mean_g_reference_ns{0};
    /// The `intersect / h_build / mwpm_build / blossom_on_h / harvest` split, as fractions of the
    /// measured total. They sum to less than 1 by exactly the unattributed cost.
    double frac_intersect{0};
    double frac_h_build{0};
    double frac_mwpm_build{0};
    double frac_blossom_on_h{0};
    double frac_harvest{0};

    double mean_degree{0};
    double mean_h_nodes{0};
    double mean_h_edges{0};
    double mean_residual_density{0};
    double mean_restarts{0};
    double q{0};
    /// §M2.6 rates, recorded in the exit artifact.
    double residual_tie_rate{0};
    double pairing_tie_rate{0};
    double mask_divergence_rate{0};
};

inline BallSummary summarize_ball(const BallAggregateStats& stats) {
    BallSummary summary;
    if (stats.shots == 0)
        return summary;
    double shots = (double)stats.shots;
    if (stats.sum_total_ns > 0) {
        summary.speedup_vs_m1 = (double)stats.sum_g_reference_ns / (double)stats.sum_total_ns;
        double total = (double)stats.sum_total_ns;
        summary.frac_intersect = (double)stats.sum_intersect_ns / total;
        summary.frac_h_build = (double)stats.sum_h_build_ns / total;
        summary.frac_mwpm_build = (double)stats.sum_mwpm_build_ns / total;
        summary.frac_blossom_on_h = (double)stats.sum_blossom_on_h_ns / total;
        summary.frac_harvest = (double)stats.sum_harvest_ns / total;
    }
    summary.mean_total_ns = (double)stats.sum_total_ns / shots;
    summary.mean_g_reference_ns = (double)stats.sum_g_reference_ns / shots;
    summary.p50_total_ns = percentile_of(stats.per_shot_total_ns, 0.5);
    summary.p99_total_ns = percentile_of(stats.per_shot_total_ns, 0.99);
    summary.p999_total_ns = percentile_of(stats.per_shot_total_ns, 0.999);
    summary.mean_h_nodes = (double)stats.sum_h_nodes / shots;
    summary.mean_h_edges = (double)stats.sum_h_edges / shots;
    if (stats.sum_h_nodes)
        summary.mean_degree = 2.0 * (double)stats.sum_h_edges / (double)stats.sum_h_nodes;
    if (stats.sum_n_defects)
        summary.mean_residual_density = (double)stats.sum_residual_size / (double)stats.sum_n_defects;
    summary.mean_restarts = (double)stats.sum_restarts / shots;
    summary.q = (double)stats.shots_truncated / shots;
    summary.residual_tie_rate = (double)stats.residual_ties / shots;
    summary.pairing_tie_rate = (double)stats.pairing_ties / shots;
    summary.mask_divergence_rate = (double)stats.mask_divergences / shots;
    return summary;
}

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_PERF_BALL_PROFILE_H
