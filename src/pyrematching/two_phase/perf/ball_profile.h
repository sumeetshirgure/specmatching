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
    /// The solve on `H`: `process_timeline_until_horizon` on the truncated path, or stock blossom
    /// run to completion on §M7's. On the latter it excludes `dual_scan_ns`, so the two are
    /// additive and the certificate's cost is visible on its own.
    long long blossom_on_h_ns{0};
    long long harvest_ns{0};
    /// §M7.7. The one terminal `max_u Y(u)` scan, paid only on shots that completed on `H`. This is
    /// the only cost the certificate scheme adds over an unmodified stock-on-`H` solve, and the
    /// §M7.8 read subtracts it from the gating overhead the swap removes.
    long long dual_scan_ns{0};
    /// §M3.4: tearing down an escalating shot's `H` instance without harvesting it. Non-zero only
    /// on shots that truncate, and included in `harvest_ns`'s window, which is why it is reported
    /// beside it rather than added to it.
    long long abandon_ns{0};
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

    /// §M7.7's per-shot certificate readout. All three are zero on the truncated-`H` path, where
    /// there is no certificate to report; on §M7's path exactly one of `certified` and
    /// "escalates" holds, and `h_no_perfect_matching` says which of the two escalation triggers
    /// fired.
    ///
    /// 1 iff the shot completed on `H` **and** `max_u Y(u) <= T_int` — i.e. iff the shot is
    /// certified globally optimal and `H`'s answer is kept.
    int certified{0};
    /// 1 iff `H` admitted no perfect matching. A valid escalation trigger, not a failure (§M7.0).
    int h_no_perfect_matching{0};
    /// `max_u Y(u)` at completion: the quantity the certificate tests. Undefined, and left at zero,
    /// when `H` could not complete.
    pm::total_weight_int max_dual_at_completion{0};

    /// §M2 structural (hardware-budget) counters. `isect_ns`, `h_build_ns` and `mwpm_build_ns`
    /// measure this laptop's DRAM latency; these measure what the three stages we intend to move
    /// off the CPU actually *move*, which is machine-independent and can be fed to a latency model
    /// for the target architecture.
    ///
    /// `isect_bytes` is build-mode dependent, hence `isect_mode_bitset` next to it, and
    /// `isect_scan_bytes_other_mode` — what the other mode's traversal would have read on this
    /// same shot, derived from the tables rather than measured by running it.
    uint64_t isect_scan_bytes{0};
    uint64_t isect_hit_bytes{0};
    uint64_t isect_scan_bytes_other_mode{0};
    bool isect_mode_bitset{false};
    /// Edge records `H` construction emits: undirected defect-defect pairs, each once, plus
    /// boundary edges. This is the crossbar's area; `max_degree` is its width.
    int hbld_edges_written{0};
    /// What `BallMwpm::rebuild` writes, counted at the write sites. See `BallMwpmCounts`.
    int mwpm_init_node_elements{0};
    int mwpm_init_edge_elements{0};

    inline uint64_t isect_bytes() const {
        return isect_scan_bytes + isect_hit_bytes;
    }
    inline int mwpm_init_elements() const {
        return mwpm_init_node_elements + mwpm_init_edge_elements;
    }

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
    int boundary_ties{0};

    /// §M2.9.6. The split of `harvest_ns` into its four stages: enumerating the live regions and
    /// alternating tree nodes, the base descent into exposed root blossoms, shattering matched
    /// blossoms, and the reductions plus residual compaction. This is the software A/B that says
    /// how much of harvest's 7–17% of the total was enumeration rather than extraction.
    long long harvest_enumerate_ns{0};
    long long harvest_reduce_ns{0};
    long long harvest_base_descent_ns{0};
    long long harvest_shatter_ns{0};

    /// §M2.9.6 measurements 2 and 3, in counts and depths rather than nanoseconds.
    ///
    /// `blossom_formations` is the number of blossoms created during this shot's solve — the
    /// quantity §M2.9.4's eager-cached-base option has to be rarer than base readouts to be worth
    /// its field and its mutation site. The other two describe the blossoms alive at truncation.
    int blossom_formations{0};
    int max_blossom_nesting_depth{0};
    int max_blossom_members{0};
    int matched_blossom_shatters{0};
    /// The exposed-root-blossom restriction of the two above. §M2.9.4's base descent walks these
    /// and only these, so these are the numbers that choose its mechanism.
    int max_exposed_blossom_depth{0};
    int max_exposed_blossom_members{0};

    /// §M2.9.6 measurements 1 and 4, and harvest's own modelled depth. Profiling builds only.
    int largest_tree_size{0};
    int harvest_dependent_depth{0};
    int solve_dependent_depth{0};
    int solve_events{0};

    /// §C.1's component structure of this shot's `H`. Filled by
    /// `BallDecoder::analyze_last_shot_components`, which runs **after** `total_ns` has been read
    /// and is timed by nothing: the component work is assumed free (hardware-offloadable), so it is
    /// never charged to `blossom_on_h_ns`, `harvest_ns`, `total_ns` or any other reported latency.
    ///
    /// Left at `measured == 0` unless `BallConfig::collect_component_stats` is on, so a campaign
    /// that did not collect it reports that rather than a row of zeros.
    ComponentStats components;

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
    /// §M7.7. Zero on the truncated-`H` path; on §M7's path `shots_certified` and
    /// `shots - shots_certified` are the kept and the escalating shots respectively.
    uint64_t shots_certified{0};
    uint64_t shots_h_no_perfect_matching{0};
    uint64_t max_dual_at_completion{0};

    long long sum_intersect_ns{0};
    long long sum_h_build_ns{0};
    long long sum_mwpm_build_ns{0};
    long long sum_blossom_on_h_ns{0};
    long long sum_harvest_ns{0};
    long long sum_dual_scan_ns{0};
    long long sum_abandon_ns{0};
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
    uint64_t boundary_ties{0};

    /// §M2 structural counters. Sums *and* maxima: the latency budget of a hardware stage is set by
    /// the worst shot it has to absorb, not by the mean.
    uint64_t sum_isect_scan_bytes{0};
    uint64_t sum_isect_hit_bytes{0};
    uint64_t sum_isect_bytes{0};
    uint64_t sum_isect_scan_bytes_other_mode{0};
    uint64_t max_isect_scan_bytes{0};
    uint64_t max_isect_hit_bytes{0};
    uint64_t max_isect_bytes{0};
    uint64_t max_isect_scan_bytes_other_mode{0};
    uint64_t sum_hbld_edges_written{0};
    uint64_t max_hbld_edges_written{0};
    uint64_t sum_mwpm_init_node_elements{0};
    uint64_t sum_mwpm_init_edge_elements{0};
    uint64_t sum_mwpm_init_elements{0};
    uint64_t max_mwpm_init_elements{0};
    /// Shots decoded in `BITSET` mode. Equal to `shots` or 0 in any single-mode campaign, which is
    /// how `BallSummary::mode_bitset` decides what to label the row.
    uint64_t shots_bitset_mode{0};

    /// §M2.9.6 accumulators.
    long long sum_harvest_enumerate_ns{0};
    long long sum_harvest_reduce_ns{0};
    long long sum_harvest_base_descent_ns{0};
    long long sum_harvest_shatter_ns{0};
    uint64_t sum_blossom_formations{0};
    uint64_t max_blossom_nesting_depth{0};
    uint64_t max_blossom_members{0};
    uint64_t max_exposed_blossom_depth{0};
    uint64_t max_exposed_blossom_members{0};
    /// Histogram of the exposed-root-blossom nesting depth, capped, last bin overflowing. This is
    /// the distribution the M2.9 exit checkpoint asks be written down next to §M2.9.4's decision.
    static constexpr size_t DEPTH_HIST_BINS = 9;
    std::vector<uint64_t> exposed_depth_hist = std::vector<uint64_t>(DEPTH_HIST_BINS, 0);
    uint64_t sum_matched_blossom_shatters{0};
    uint64_t max_largest_tree_size{0};
    uint64_t sum_largest_tree_size{0};
    uint64_t max_harvest_dependent_depth{0};
    uint64_t sum_harvest_dependent_depth{0};
    uint64_t max_solve_dependent_depth{0};
    uint64_t sum_solve_dependent_depth{0};
    uint64_t sum_solve_events{0};

    bool keep_per_shot{false};
    std::vector<long long> per_shot_total_ns;
    std::vector<long long> per_shot_g_reference_ns;
    /// Benchmark mode only, exactly like the wall times above, so the structural counters get the
    /// same p99/p999 treatment.
    std::vector<long long> per_shot_isect_bytes;
    std::vector<long long> per_shot_hbld_edges_written;
    std::vector<long long> per_shot_mwpm_init_elements;

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
        sum_dual_scan_ns += profile.dual_scan_ns;
        sum_abandon_ns += profile.abandon_ns;
        shots_certified += (uint64_t)profile.certified;
        shots_h_no_perfect_matching += (uint64_t)profile.h_no_perfect_matching;
        max_dual_at_completion = std::max(max_dual_at_completion, (uint64_t)profile.max_dual_at_completion);
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
        boundary_ties += (uint64_t)profile.boundary_ties;

        sum_isect_scan_bytes += profile.isect_scan_bytes;
        sum_isect_hit_bytes += profile.isect_hit_bytes;
        sum_isect_bytes += profile.isect_bytes();
        sum_isect_scan_bytes_other_mode += profile.isect_scan_bytes_other_mode;
        max_isect_scan_bytes = std::max(max_isect_scan_bytes, profile.isect_scan_bytes);
        max_isect_hit_bytes = std::max(max_isect_hit_bytes, profile.isect_hit_bytes);
        max_isect_bytes = std::max(max_isect_bytes, profile.isect_bytes());
        max_isect_scan_bytes_other_mode =
            std::max(max_isect_scan_bytes_other_mode, profile.isect_scan_bytes_other_mode);
        sum_hbld_edges_written += (uint64_t)profile.hbld_edges_written;
        max_hbld_edges_written = std::max(max_hbld_edges_written, (uint64_t)profile.hbld_edges_written);
        sum_mwpm_init_node_elements += (uint64_t)profile.mwpm_init_node_elements;
        sum_mwpm_init_edge_elements += (uint64_t)profile.mwpm_init_edge_elements;
        sum_mwpm_init_elements += (uint64_t)profile.mwpm_init_elements();
        max_mwpm_init_elements = std::max(max_mwpm_init_elements, (uint64_t)profile.mwpm_init_elements());
        shots_bitset_mode += profile.isect_mode_bitset ? 1 : 0;

        sum_harvest_enumerate_ns += profile.harvest_enumerate_ns;
        sum_harvest_reduce_ns += profile.harvest_reduce_ns;
        sum_harvest_base_descent_ns += profile.harvest_base_descent_ns;
        sum_harvest_shatter_ns += profile.harvest_shatter_ns;
        sum_blossom_formations += (uint64_t)profile.blossom_formations;
        max_blossom_nesting_depth = std::max(max_blossom_nesting_depth, (uint64_t)profile.max_blossom_nesting_depth);
        max_blossom_members = std::max(max_blossom_members, (uint64_t)profile.max_blossom_members);
        max_exposed_blossom_depth = std::max(max_exposed_blossom_depth, (uint64_t)profile.max_exposed_blossom_depth);
        max_exposed_blossom_members =
            std::max(max_exposed_blossom_members, (uint64_t)profile.max_exposed_blossom_members);
        if (profile.max_exposed_blossom_depth > 0) {
            size_t bin = std::min((size_t)profile.max_exposed_blossom_depth, DEPTH_HIST_BINS - 1);
            exposed_depth_hist[bin]++;
        }
        sum_matched_blossom_shatters += (uint64_t)profile.matched_blossom_shatters;
        max_largest_tree_size = std::max(max_largest_tree_size, (uint64_t)profile.largest_tree_size);
        sum_largest_tree_size += (uint64_t)profile.largest_tree_size;
        max_harvest_dependent_depth = std::max(max_harvest_dependent_depth, (uint64_t)profile.harvest_dependent_depth);
        sum_harvest_dependent_depth += (uint64_t)profile.harvest_dependent_depth;
        max_solve_dependent_depth = std::max(max_solve_dependent_depth, (uint64_t)profile.solve_dependent_depth);
        sum_solve_dependent_depth += (uint64_t)profile.solve_dependent_depth;
        sum_solve_events += (uint64_t)profile.solve_events;

        if (keep_per_shot) {
            per_shot_total_ns.push_back(profile.total_ns);
            per_shot_g_reference_ns.push_back(profile.g_reference_ns);
            per_shot_isect_bytes.push_back((long long)profile.isect_bytes());
            per_shot_hbld_edges_written.push_back((long long)profile.hbld_edges_written);
            per_shot_mwpm_init_elements.push_back((long long)profile.mwpm_init_elements());
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
    /// §M7.7/§M7.8: the certificate's own share of the shot, so the stage split of the M7 read
    /// (`intersect / h_build / mwpm_build / blossom_on_h / dual_scan`) is complete.
    double frac_dual_scan{0};
    double mean_dual_scan_ns{0};
    /// §M7: fraction of shots kept on `H` under the certificate, and the share of the escalating
    /// ones that escalated because `H` had no perfect matching rather than because the terminal
    /// dual escaped `T`. Both zero on the truncated-`H` path.
    double certified_rate{0};
    double h_no_perfect_matching_rate{0};

    double mean_degree{0};
    double max_degree{0};
    double mean_h_nodes{0};
    double mean_h_edges{0};
    double mean_residual_density{0};
    double mean_restarts{0};
    double q{0};
    /// §M2.6 rates, recorded in the exit artifact.
    double residual_tie_rate{0};
    double pairing_tie_rate{0};
    double boundary_tie_rate{0};
    double mask_divergence_rate{0};

    /// §M2.9.6. Harvest's four stages as fractions of measured `harvest_ns`. The A/B the design
    /// asks for is `frac_harvest_enumerate`: if most of harvest was enumeration rather than
    /// extraction, §M2.9.1–§M2.9.2 are a CPU win as well as a critical-path one.
    double frac_harvest_enumerate{0};
    double frac_harvest_reduce{0};
    double frac_harvest_base_descent{0};
    double frac_harvest_shatter{0};

    /// Measurements 1–4, in counts and depths.
    double mean_blossom_formations{0};
    double mean_matched_blossom_shatters{0};
    double mean_largest_tree_size{0};
    double max_largest_tree_size{0};
    double max_blossom_nesting_depth{0};
    double mean_harvest_dependent_depth{0};
    double max_harvest_dependent_depth{0};
    double mean_solve_dependent_depth{0};
    double max_solve_dependent_depth{0};
    double mean_solve_events{0};
    /// The number the M2.9 exit checkpoint asks for: harvest's share of the critical path that is
    /// left once the embarrassingly parallel stages collapse, measured rather than estimated.
    double harvest_share_of_critical_path{0};

    /// §M2 structural counters, per shot unless the name says otherwise.
    ///
    /// `isect_bytes_per_defect` is the one that gets quoted: it is the local memory a per-defect
    /// processing element has to be able to stream, which is the M2 exit checkpoint's headline
    /// number for the local-memory hardware argument. The divisor is `h_nodes` — the post-preamble
    /// defects, which is exactly the loop's trip count — not the raw detection-event count.
    bool mode_bitset{false};
    double mean_isect_bytes{0};
    double mean_isect_scan_bytes{0};
    double mean_isect_hit_bytes{0};
    double mean_isect_bytes_other_mode{0};
    double isect_bytes_per_defect{0};
    double isect_scan_bytes_per_defect{0};
    double isect_hit_bytes_per_defect{0};
    double isect_bytes_other_mode_per_defect{0};
    double max_isect_bytes{0};
    double max_isect_scan_bytes{0};
    double max_isect_hit_bytes{0};
    double max_isect_bytes_other_mode{0};
    double p99_isect_bytes{0};
    double p999_isect_bytes{0};
    double mean_hbld_edges_written{0};
    double max_hbld_edges_written{0};
    double p99_hbld_edges_written{0};
    double p999_hbld_edges_written{0};
    double mean_mwpm_init_elements{0};
    double mean_mwpm_init_node_elements{0};
    double mean_mwpm_init_edge_elements{0};
    double max_mwpm_init_elements{0};
    double p99_mwpm_init_elements{0};
    double p999_mwpm_init_elements{0};
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
        summary.frac_dual_scan = (double)stats.sum_dual_scan_ns / total;
    }
    summary.mean_dual_scan_ns = (double)stats.sum_dual_scan_ns / shots;
    summary.certified_rate = (double)stats.shots_certified / shots;
    summary.h_no_perfect_matching_rate = (double)stats.shots_h_no_perfect_matching / shots;
    summary.mean_total_ns = (double)stats.sum_total_ns / shots;
    summary.mean_g_reference_ns = (double)stats.sum_g_reference_ns / shots;
    summary.p50_total_ns = percentile_of(stats.per_shot_total_ns, 0.5);
    summary.p99_total_ns = percentile_of(stats.per_shot_total_ns, 0.99);
    summary.p999_total_ns = percentile_of(stats.per_shot_total_ns, 0.999);
    summary.mean_h_nodes = (double)stats.sum_h_nodes / shots;
    summary.mean_h_edges = (double)stats.sum_h_edges / shots;
    summary.max_degree = (double)stats.max_degree;
    if (stats.sum_h_nodes)
        summary.mean_degree = 2.0 * (double)stats.sum_h_edges / (double)stats.sum_h_nodes;
    if (stats.sum_n_defects)
        summary.mean_residual_density = (double)stats.sum_residual_size / (double)stats.sum_n_defects;
    summary.mean_restarts = (double)stats.sum_restarts / shots;
    summary.q = (double)stats.shots_truncated / shots;
    summary.residual_tie_rate = (double)stats.residual_ties / shots;
    summary.pairing_tie_rate = (double)stats.pairing_ties / shots;
    summary.boundary_tie_rate = (double)stats.boundary_ties / shots;
    summary.mask_divergence_rate = (double)stats.mask_divergences / shots;

    if (stats.sum_harvest_ns > 0) {
        double harvest = (double)stats.sum_harvest_ns;
        summary.frac_harvest_enumerate = (double)stats.sum_harvest_enumerate_ns / harvest;
        summary.frac_harvest_reduce = (double)stats.sum_harvest_reduce_ns / harvest;
        summary.frac_harvest_base_descent = (double)stats.sum_harvest_base_descent_ns / harvest;
        summary.frac_harvest_shatter = (double)stats.sum_harvest_shatter_ns / harvest;
    }
    summary.mean_blossom_formations = (double)stats.sum_blossom_formations / shots;
    summary.mean_matched_blossom_shatters = (double)stats.sum_matched_blossom_shatters / shots;
    summary.mean_largest_tree_size = (double)stats.sum_largest_tree_size / shots;
    summary.max_largest_tree_size = (double)stats.max_largest_tree_size;
    summary.max_blossom_nesting_depth = (double)stats.max_blossom_nesting_depth;
    summary.mean_harvest_dependent_depth = (double)stats.sum_harvest_dependent_depth / shots;
    summary.max_harvest_dependent_depth = (double)stats.max_harvest_dependent_depth;
    summary.mean_solve_dependent_depth = (double)stats.sum_solve_dependent_depth / shots;
    summary.max_solve_dependent_depth = (double)stats.max_solve_dependent_depth;
    summary.mean_solve_events = (double)stats.sum_solve_events / shots;
    double critical_path = summary.mean_harvest_dependent_depth + summary.mean_solve_dependent_depth;
    if (critical_path > 0)
        summary.harvest_share_of_critical_path = summary.mean_harvest_dependent_depth / critical_path;

    summary.mode_bitset = stats.shots_bitset_mode == stats.shots;
    summary.mean_isect_bytes = (double)stats.sum_isect_bytes / shots;
    summary.mean_isect_scan_bytes = (double)stats.sum_isect_scan_bytes / shots;
    summary.mean_isect_hit_bytes = (double)stats.sum_isect_hit_bytes / shots;
    summary.mean_isect_bytes_other_mode = (double)stats.sum_isect_scan_bytes_other_mode / shots;
    if (stats.sum_h_nodes) {
        double defects = (double)stats.sum_h_nodes;
        summary.isect_bytes_per_defect = (double)stats.sum_isect_bytes / defects;
        summary.isect_scan_bytes_per_defect = (double)stats.sum_isect_scan_bytes / defects;
        summary.isect_hit_bytes_per_defect = (double)stats.sum_isect_hit_bytes / defects;
        summary.isect_bytes_other_mode_per_defect = (double)stats.sum_isect_scan_bytes_other_mode / defects;
    }
    summary.max_isect_bytes = (double)stats.max_isect_bytes;
    summary.max_isect_scan_bytes = (double)stats.max_isect_scan_bytes;
    summary.max_isect_hit_bytes = (double)stats.max_isect_hit_bytes;
    summary.max_isect_bytes_other_mode = (double)stats.max_isect_scan_bytes_other_mode;
    summary.p99_isect_bytes = percentile_of(stats.per_shot_isect_bytes, 0.99);
    summary.p999_isect_bytes = percentile_of(stats.per_shot_isect_bytes, 0.999);
    summary.mean_hbld_edges_written = (double)stats.sum_hbld_edges_written / shots;
    summary.max_hbld_edges_written = (double)stats.max_hbld_edges_written;
    summary.p99_hbld_edges_written = percentile_of(stats.per_shot_hbld_edges_written, 0.99);
    summary.p999_hbld_edges_written = percentile_of(stats.per_shot_hbld_edges_written, 0.999);
    summary.mean_mwpm_init_elements = (double)stats.sum_mwpm_init_elements / shots;
    summary.mean_mwpm_init_node_elements = (double)stats.sum_mwpm_init_node_elements / shots;
    summary.mean_mwpm_init_edge_elements = (double)stats.sum_mwpm_init_edge_elements / shots;
    summary.max_mwpm_init_elements = (double)stats.max_mwpm_init_elements;
    summary.p99_mwpm_init_elements = percentile_of(stats.per_shot_mwpm_init_elements, 0.99);
    summary.p999_mwpm_init_elements = percentile_of(stats.per_shot_mwpm_init_elements, 0.999);
    return summary;
}

}  // namespace two_phase
}  // namespace pm

#endif  // PYREMATCHING_TWO_PHASE_PERF_BALL_PROFILE_H
