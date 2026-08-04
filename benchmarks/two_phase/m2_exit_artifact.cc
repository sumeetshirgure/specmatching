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

/// The M2 exit artifact.
///
/// Three sections, matching the M2 exit checkpoint:
///
///  - `ball`   — ball table cost per `d`: bytes broken down per array, mean/max ball size, mean
///               bitset window, compile wall time. This is the headline number for the
///               local-memory hardware argument.
///  - `profile`— `speedup_vs_m1` per `(d, p, T)` with the
///               `intersect / h_build / mwpm_build / blossom_on_h / harvest` split, `H`'s size and
///               degree, and the §M2.6 tie and mask-divergence rates. **Expect the win to grow as
///               `p` falls**, since the `beta * V` floor that M2 removes dominates at low defect
///               density — benchmarking only at `p = 5e-3` understates the result.
///               It also carries `mean_stock_ns`: stock sparse blossom run to *completion* on `G`,
///               the un-truncated baseline. That is a different job from either truncated number —
///               it returns a full matching, where Phase 1 returns a partial one plus a residual
///               that still has to be cleaned up — so it is a separate column and never enters
///               `speedup_vs_m1`.
///  - `ler`    — level 3: the full pipeline (front end + exact cleanup of the residual) for both
///               front ends, against stock exact decoding, so that `d_eff/d` can be fitted per
///               M6.4 and compared between M1 and M2.
///
/// Horizons are quoted in multiples of one lattice edge weight, the normalisation under which the
/// M1 results are `d`-independent.
///
/// Usage:
///   two_phase_m2_artifact [--distances 5,7,9,11,13] [--error-rates 0.0005,0.001,0.003,0.005]
///                         [--horizons 1.5,2.0] [--shots 2000] [--ler-shots 20000]
///                         [--mode scan|bitset] [--seed N] [--csv path]

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/sparse_blossom/driver/user_graph.h"
#include "pyrematching/two_phase/manifold/ball_decoding.h"
#include "pyrematching/two_phase/perf/ball_profile.h"
#include "stim.h"

using namespace pm::two_phase;

namespace {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 10001;

struct Options {
    std::vector<size_t> distances = {5, 7, 9, 11, 13};
    std::vector<double> error_rates = {0.0005, 0.001, 0.003, 0.005};
    /// In multiples of one lattice edge weight.
    std::vector<double> horizons = {1.5, 2.0};
    size_t shots = 2000;
    size_t ler_shots = 20000;
    BallGraphBuildMode mode = BallGraphBuildMode::SCAN;
    uint64_t seed = 20260803;
    std::string csv_path;
};

template <typename T>
std::vector<T> parse_list(const std::string& text) {
    std::vector<T> out;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty())
            continue;
        std::stringstream item_stream(item);
        T value;
        item_stream >> value;
        out.push_back(value);
    }
    return out;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; i++) {
        std::string flag = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for " + flag);
            return argv[++i];
        };
        if (flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--horizons") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--ler-shots") {
            options.ler_shots = std::stoul(next());
        } else if (flag == "--mode") {
            std::string mode = next();
            options.mode = mode == "bitset" ? BallGraphBuildMode::BITSET : BallGraphBuildMode::SCAN;
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--csv") {
            options.csv_path = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    return options;
}

struct Experiment {
    stim::DetectorErrorModel dem;
    std::vector<std::vector<uint64_t>> shots;
    /// The true logical flips per shot, needed for the level-3 logical error rate.
    std::vector<uint8_t> observable_flips;
};

Experiment generate(size_t distance, double noise, size_t num_shots, uint64_t seed) {
    stim::CircuitGenParameters gen(distance, distance, "rotated_memory_x");
    gen.after_clifford_depolarization = noise;
    gen.after_reset_flip_probability = noise;
    gen.before_measure_flip_probability = noise;
    stim::Circuit circuit = stim::generate_surface_code_circuit(gen).circuit;

    std::mt19937_64 rng(seed);
    size_t num_detectors = circuit.count_detectors();
    auto dets_obs = stim::sample_batch_detection_events<stim::MAX_BITWORD_WIDTH>(circuit, num_shots, rng);
    auto& dets = dets_obs.first;
    auto& obs = dets_obs.second;

    Experiment experiment;
    experiment.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    for (size_t shot = 0; shot < num_shots; shot++) {
        experiment.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if (dets[d][shot])
                experiment.shots.back().push_back(d);
        }
        experiment.observable_flips.push_back(obs[0][shot] ? 1 : 0);
    }
    return experiment;
}

double edge_weight_units(const pm::MatchingGraph& graph) {
    std::vector<pm::weight_int> weights;
    for (const auto& node : graph.nodes) {
        for (auto weight : node.neighbor_weights)
            weights.push_back(weight);
    }
    if (weights.empty())
        return 1.0;
    std::sort(weights.begin(), weights.end());
    return (double)weights[weights.size() / 2] / graph.normalising_constant;
}

pm::Mwpm build_mwpm(const stim::DetectorErrorModel& dem) {
    return pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
}

/// Rows are emitted to stdout as a table and, optionally, to a CSV. Each row carries its section so
/// one file holds the whole artifact.
struct CsvWriter {
    std::ofstream stream;
    bool enabled{false};

    void open(const std::string& path) {
        if (path.empty())
            return;
        stream.open(path);
        enabled = stream.is_open();
        if (!enabled)
            std::cerr << "warning: could not open " << path << " for writing\n";
    }
    void row(const std::string& text) {
        if (enabled)
            stream << text << "\n";
    }
};

/// A one-standard-error linear fit of `ln(LER)` against `d` at fixed `p`. `c2 = d_eff/d` is
/// proportional to the slope; what M2 has to reproduce is M1's slope, not an absolute value.
struct Fit {
    double slope{0};
    double stderr_slope{0};
    size_t points{0};
};

Fit fit_log_ler(const std::vector<double>& distances, const std::vector<double>& lers) {
    Fit fit;
    std::vector<double> xs;
    std::vector<double> ys;
    for (size_t i = 0; i < lers.size(); i++) {
        if (lers[i] > 0) {
            xs.push_back(distances[i]);
            ys.push_back(std::log(lers[i]));
        }
    }
    fit.points = xs.size();
    if (xs.size() < 3)
        return fit;
    double n = (double)xs.size();
    double mean_x = 0;
    double mean_y = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        mean_x += xs[i] / n;
        mean_y += ys[i] / n;
    }
    double sxx = 0;
    double sxy = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        sxx += (xs[i] - mean_x) * (xs[i] - mean_x);
        sxy += (xs[i] - mean_x) * (ys[i] - mean_y);
    }
    if (sxx == 0)
        return fit;
    fit.slope = sxy / sxx;
    double residual = 0;
    for (size_t i = 0; i < xs.size(); i++) {
        double predicted = mean_y + fit.slope * (xs[i] - mean_x);
        residual += (ys[i] - predicted) * (ys[i] - predicted);
    }
    if (xs.size() > 2)
        fit.stderr_slope = std::sqrt(residual / (n - 2) / sxx);
    return fit;
}

}  // namespace

int main(int argc, char** argv) {
    Options options;
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }

    CsvWriter csv;
    csv.open(options.csv_path);
    csv.row("section,d,p,T,key,value");
    auto emit = [&](const std::string& section, size_t d, double p, double t, const std::string& key, double value) {
        std::ostringstream row;
        row << section << "," << d << "," << p << "," << t << "," << key << "," << value;
        csv.row(row.str());
    };

    // ---------------------------------------------------------------- section: ball table cost
    std::printf("\n=== ball tables (T_max = 2, R = 4 edge weights) ===\n");
    std::printf(
        "%4s %10s %12s %12s %10s %10s %10s %10s\n",
        "d",
        "nodes",
        "entries",
        "bytes",
        "mean|B|",
        "max|B|",
        "words",
        "compile_s");
    for (size_t distance : options.distances) {
        Experiment experiment = generate(distance, options.error_rates.back(), 1, options.seed);
        pm::Mwpm mwpm = build_mwpm(experiment.dem);
        double unit = edge_weight_units(mwpm.flooder.graph);
        BallParams params;
        params.T_max = 2.0 * unit;
        params.R = 4.0 * unit;
        BallTables tables = compile_ball_tables(mwpm, params);
        const BallStats& stats = tables.stats;
        std::printf(
            "%4zu %10llu %12llu %12llu %10.1f %10llu %10.1f %10.2f\n",
            distance,
            (unsigned long long)stats.num_nodes,
            (unsigned long long)stats.total_entries,
            (unsigned long long)stats.bytes_total,
            stats.mean_ball_size,
            (unsigned long long)stats.max_ball_size,
            stats.mean_ball_word_len,
            stats.compile_wall_seconds);
        emit("ball", distance, 0, 0, "nodes", (double)stats.num_nodes);
        emit("ball", distance, 0, 0, "entries", (double)stats.total_entries);
        emit("ball", distance, 0, 0, "bytes_total", (double)stats.bytes_total);
        emit("ball", distance, 0, 0, "bytes_targets", (double)stats.bytes_targets);
        emit("ball", distance, 0, 0, "bytes_weights", (double)stats.bytes_weights);
        emit("ball", distance, 0, 0, "bytes_masks", (double)stats.bytes_masks);
        emit("ball", distance, 0, 0, "bytes_words", (double)stats.bytes_words);
        emit("ball", distance, 0, 0, "bytes_offsets", (double)stats.bytes_offsets);
        emit("ball", distance, 0, 0, "mean_ball_size", stats.mean_ball_size);
        emit("ball", distance, 0, 0, "max_ball_size", (double)stats.max_ball_size);
        emit("ball", distance, 0, 0, "mean_ball_word_len", stats.mean_ball_word_len);
        emit("ball", distance, 0, 0, "ambiguous_mask_pairs", (double)stats.ambiguous_mask_pairs);
        emit("ball", distance, 0, 0, "compile_wall_seconds", stats.compile_wall_seconds);
    }

    // ------------------------------------------------------------------- section: speedup vs M1
    //
    // The structural counters share this section's sweep — they are per-shot quantities of the same
    // decodes — but they are a different kind of number (bytes and elements, not nanoseconds), so
    // they get their own table, printed after this one. Rows are buffered rather than printed
    // inline so that the two tables do not interleave.
    std::vector<std::string> structural_rows;
    std::printf("\n=== speedup vs M1 Phase 1 on G (%zu shots per point) ===\n", options.shots);
    std::printf(
        "%4s %8s %5s %9s %9s %9s %8s %8s %9s %8s %8s %7s %7s %7s %7s %7s %8s %8s\n",
        "d",
        "p",
        "T",
        "m2_ns",
        "m1_ns",
        "stock_ns",
        "speedup",
        "spd_stk",
        "h_nodes",
        "h_edges",
        "deg",
        "isect",
        "hbuild",
        "mwpm",
        "blossom",
        "harvest",
        "res_tie",
        "pair_tie");
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, noise, options.shots, options.seed + distance * 131);
            pm::Mwpm probe = build_mwpm(experiment.dem);
            double unit = edge_weight_units(probe.flooder.graph);

            // The un-truncated baseline: stock sparse blossom decoded to completion on `G`, timed
            // over the same shots. It does not depend on `T`, so it is measured once per `(d, p)`
            // and repeated across the horizon rows. Its own pass, not interleaved with the M2 path,
            // so neither one is timed with the other's working set in cache.
            pm::Mwpm stock = build_mwpm(experiment.dem);
            long long stock_total_ns = 0;
            if (!experiment.shots.empty())
                pm::decode_detection_events_for_up_to_64_observables(stock, experiment.shots[0], false);
            for (const auto& shot : experiment.shots) {
                HiResTimer stock_timer;
                stock_timer.start();
                pm::decode_detection_events_for_up_to_64_observables(stock, shot, false);
                stock_total_ns += stock_timer.elapsed_ns();
            }
            double mean_stock_ns = (double)stock_total_ns / (double)std::max<size_t>(1, experiment.shots.size());

            for (double t_edges : options.horizons) {
                BallConfig config;
                config.T = t_edges * unit;
                config.ball.T_max = 2.0 * unit;
                config.ball.R = 4.0 * unit;
                config.mode = options.mode;
                // The oracle path is the *point* of this section: `g_reference_ns` is M1's Phase 1
                // on the same shot, and it also collects the §M2.6 tie rates.
                config.verify_against_g = true;
                BallDecoder decoder = BallDecoder::from_mwpm(build_mwpm(experiment.dem), config);

                BallAggregateStats stats;
                stats.keep_per_shot = true;
                BallProfile profile;
                // One warmup shot, so the arena and the H-side Mwpm are at steady state and the
                // timings are not dominated by their first growth.
                if (!experiment.shots.empty())
                    decoder.decode_phase1(experiment.shots[0], &profile);
                for (const auto& shot : experiment.shots) {
                    HarvestResult harvest = decoder.decode_phase1(shot, &profile);
                    stats.accumulate(profile, harvest);
                }
                BallSummary summary = summarize_ball(stats);

                // The §M2.6 tie rates, collected in a separate *untimed* pass. The committed pair
                // set only exists in the match-edges flavour, so the pairing rate cannot be read
                // off the obs-flavour hot path above without changing what is being timed.
                //
                // The §M2 structural counters ride along on this pass for the same reason, one
                // step stronger: charging every pair its observable id list costs 26-44% of
                // `intersect_ns`, so collecting them in the timed pass would corrupt the very
                // split printed above. They are structural — a function of the shot and the
                // tables — so measuring them here rather than there changes nothing about them.
                decoder.config.collect_structural_counters = true;
                BallAggregateStats structural_stats;
                structural_stats.keep_per_shot = true;
                BallProfile tie_profile;
                std::vector<CommittedPair> pairs;
                uint64_t residual_ties = 0;
                uint64_t pairing_ties = 0;
                for (const auto& shot : experiment.shots) {
                    HarvestResult harvest = decoder.decode_phase1_to_match_edges(shot, pairs, &tie_profile);
                    structural_stats.accumulate(tie_profile, harvest);
                    residual_ties += (uint64_t)tie_profile.residual_ties;
                    pairing_ties += (uint64_t)tie_profile.pairing_ties;
                }
                decoder.config.collect_structural_counters = false;
                BallSummary structural = summarize_ball(structural_stats);
                double tie_shots = (double)std::max<size_t>(1, experiment.shots.size());
                summary.residual_tie_rate = (double)residual_ties / tie_shots;
                summary.pairing_tie_rate = (double)pairing_ties / tie_shots;

                std::printf(
                    "%4zu %8.4f %5.2f %9.0f %9.0f %9.0f %8.2f %8.2f %9.1f %8.1f %8.2f %7.3f %7.3f %7.3f %7.3f %7.3f "
                    "%8.4f %8.4f\n",
                    distance,
                    noise,
                    t_edges,
                    summary.mean_total_ns,
                    summary.mean_g_reference_ns,
                    mean_stock_ns,
                    summary.speedup_vs_m1,
                    summary.mean_total_ns > 0 ? mean_stock_ns / summary.mean_total_ns : 0.0,
                    summary.mean_h_nodes,
                    summary.mean_h_edges,
                    summary.mean_degree,
                    summary.frac_intersect,
                    summary.frac_h_build,
                    summary.frac_mwpm_build,
                    summary.frac_blossom_on_h,
                    summary.frac_harvest,
                    summary.residual_tie_rate,
                    summary.pairing_tie_rate);

                emit("profile", distance, noise, t_edges, "mean_total_ns", summary.mean_total_ns);
                emit("profile", distance, noise, t_edges, "mean_g_reference_ns", summary.mean_g_reference_ns);
                emit("profile", distance, noise, t_edges, "speedup_vs_m1", summary.speedup_vs_m1);
                emit("profile", distance, noise, t_edges, "mean_stock_ns", mean_stock_ns);
                emit("profile", distance, noise, t_edges, "p50_total_ns", summary.p50_total_ns);
                emit("profile", distance, noise, t_edges, "p99_total_ns", summary.p99_total_ns);
                emit("profile", distance, noise, t_edges, "p999_total_ns", summary.p999_total_ns);
                emit("profile", distance, noise, t_edges, "mean_h_nodes", summary.mean_h_nodes);
                emit("profile", distance, noise, t_edges, "mean_h_edges", summary.mean_h_edges);
                emit("profile", distance, noise, t_edges, "mean_degree", summary.mean_degree);
                emit("profile", distance, noise, t_edges, "max_degree", (double)stats.max_degree);
                emit("profile", distance, noise, t_edges, "frac_intersect", summary.frac_intersect);
                emit("profile", distance, noise, t_edges, "frac_h_build", summary.frac_h_build);
                emit("profile", distance, noise, t_edges, "frac_mwpm_build", summary.frac_mwpm_build);
                emit("profile", distance, noise, t_edges, "frac_blossom_on_h", summary.frac_blossom_on_h);
                emit("profile", distance, noise, t_edges, "frac_harvest", summary.frac_harvest);
                emit("profile", distance, noise, t_edges, "residual_density", summary.mean_residual_density);
                emit("profile", distance, noise, t_edges, "q", summary.q);
                emit("profile", distance, noise, t_edges, "residual_tie_rate", summary.residual_tie_rate);
                emit("profile", distance, noise, t_edges, "pairing_tie_rate", summary.pairing_tie_rate);
                emit("profile", distance, noise, t_edges, "mean_restarts", summary.mean_restarts);

                // The §M2 structural counters: what the three stages we intend to move off the CPU
                // — intersect, H build, MWPM rebuild — actually move, in bytes and element writes
                // rather than in this laptop's nanoseconds.
                char row[512];
                std::snprintf(
                    row,
                    sizeof(row),
                    "%4zu %8.4f %5.2f %7s %9.1f %8.1f %12.1f %10.1f %12.1f %10.1f %12.1f %9.2f %6.0f %10.2f",
                    distance,
                    noise,
                    t_edges,
                    structural.mode_bitset ? "bitset" : "scan",
                    structural.mean_h_nodes,
                    structural.mean_h_edges,
                    structural.mean_isect_bytes,
                    structural.isect_bytes_per_defect,
                    structural.mean_isect_scan_bytes,
                    structural.mean_isect_hit_bytes,
                    structural.mean_isect_bytes_other_mode,
                    structural.mean_hbld_edges_written,
                    structural.max_degree,
                    structural.mean_mwpm_init_elements);
                structural_rows.emplace_back(row);

                emit("profile", distance, noise, t_edges, "isect_mode_bitset", structural.mode_bitset ? 1.0 : 0.0);
                emit("profile", distance, noise, t_edges, "isect_bytes", structural.mean_isect_bytes);
                emit("profile", distance, noise, t_edges, "isect_scan_bytes", structural.mean_isect_scan_bytes);
                emit("profile", distance, noise, t_edges, "isect_hit_bytes", structural.mean_isect_hit_bytes);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "isect_bytes_other_mode",
                    structural.mean_isect_bytes_other_mode);
                emit("profile", distance, noise, t_edges, "isect_bytes_per_defect", structural.isect_bytes_per_defect);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "isect_scan_bytes_per_defect",
                    structural.isect_scan_bytes_per_defect);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "isect_hit_bytes_per_defect",
                    structural.isect_hit_bytes_per_defect);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "isect_bytes_other_mode_per_defect",
                    structural.isect_bytes_other_mode_per_defect);
                emit("profile", distance, noise, t_edges, "max_isect_bytes", structural.max_isect_bytes);
                emit("profile", distance, noise, t_edges, "max_isect_scan_bytes", structural.max_isect_scan_bytes);
                emit("profile", distance, noise, t_edges, "max_isect_hit_bytes", structural.max_isect_hit_bytes);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "max_isect_bytes_other_mode",
                    structural.max_isect_bytes_other_mode);
                emit("profile", distance, noise, t_edges, "p99_isect_bytes", structural.p99_isect_bytes);
                emit("profile", distance, noise, t_edges, "p999_isect_bytes", structural.p999_isect_bytes);
                emit("profile", distance, noise, t_edges, "hbld_edges_written", structural.mean_hbld_edges_written);
                emit("profile", distance, noise, t_edges, "max_hbld_edges_written", structural.max_hbld_edges_written);
                emit("profile", distance, noise, t_edges, "p99_hbld_edges_written", structural.p99_hbld_edges_written);
                emit(
                    "profile", distance, noise, t_edges, "p999_hbld_edges_written", structural.p999_hbld_edges_written);
                emit("profile", distance, noise, t_edges, "h_max_degree", structural.max_degree);
                emit("profile", distance, noise, t_edges, "mwpm_init_elements", structural.mean_mwpm_init_elements);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "mwpm_init_node_elements",
                    structural.mean_mwpm_init_node_elements);
                emit(
                    "profile",
                    distance,
                    noise,
                    t_edges,
                    "mwpm_init_edge_elements",
                    structural.mean_mwpm_init_edge_elements);
                emit("profile", distance, noise, t_edges, "max_mwpm_init_elements", structural.max_mwpm_init_elements);
                emit("profile", distance, noise, t_edges, "p99_mwpm_init_elements", structural.p99_mwpm_init_elements);
                emit(
                    "profile", distance, noise, t_edges, "p999_mwpm_init_elements", structural.p999_mwpm_init_elements);
            }
        }
    }

    // -------------------------------------------------------- section: structural (budget) counters
    //
    // Wall time for `isect`, `hbld` and `mwpm` measures this laptop's DRAM latency, not the target
    // architecture. These are the machine-independent quantities behind those three stages: bytes
    // moved out of the ball tables, edge records emitted, element writes into the `Mwpm` on `H`.
    // `isect_B/def` is the per-PE local-memory budget; `other_B` is what the *other* build mode's
    // traversal would have read on the same shots, so the SCAN/BITSET crossover can be read off
    // without running both.
    std::printf("\n=== structural counters, per shot (%zu shots per point) ===\n", options.shots);
    std::printf(
        "%4s %8s %5s %7s %9s %8s %12s %10s %12s %10s %12s %9s %6s %10s\n",
        "d",
        "p",
        "T",
        "mode",
        "h_nodes",
        "h_edges",
        "isect_B",
        "isect_B/def",
        "scan_B",
        "hit_B",
        "other_B",
        "hbld_edges",
        "maxdeg",
        "mwpm_init");
    for (const std::string& row : structural_rows)
        std::printf("%s\n", row.c_str());
    std::printf(
        "element sizes (bytes): ball_target %llu, ball_w_int %llu, ball_words %llu, ball_word_rank %llu, "
        "ball_entry_by_rank %llu, ball_mask_offsets %llu, ball_mask_ids %llu, has_bcost %llu, bcost_w_int %llu\n",
        (unsigned long long)ball_element_bytes::TARGET,
        (unsigned long long)ball_element_bytes::WEIGHT,
        (unsigned long long)ball_element_bytes::WORD,
        (unsigned long long)ball_element_bytes::WORD_RANK,
        (unsigned long long)ball_element_bytes::ENTRY_BY_RANK,
        (unsigned long long)ball_element_bytes::MASK_OFFSET,
        (unsigned long long)ball_element_bytes::MASK_ID,
        (unsigned long long)ball_element_bytes::HAS_BCOST,
        (unsigned long long)ball_element_bytes::BCOST_WEIGHT);

    // ------------------------------------------------- section: level 3, end-to-end logical error
    //
    // Full pipeline for both front ends: Phase 1 (truncated blossom, on `G` for M1 and on `H` for
    // M2) followed by an *exact* decode of the residual on `G`. Stock exact decoding of the whole
    // syndrome is the third curve. `d_eff/d` is proportional to the slope of `ln(LER)` against `d`,
    // and what M2 must reproduce is M1's slope.
    std::printf("\n=== level 3: end-to-end LER at T = 2 (%zu shots per point) ===\n", options.ler_shots);
    std::printf("%4s %8s %12s %12s %12s\n", "d", "p", "ler_exact", "ler_m1", "ler_m2");
    std::map<double, std::vector<double>> ler_exact_by_p;
    std::map<double, std::vector<double>> ler_m1_by_p;
    std::map<double, std::vector<double>> ler_m2_by_p;
    std::map<double, std::vector<double>> distances_by_p;

    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, noise, options.ler_shots, options.seed + 977 * distance);
            pm::Mwpm probe = build_mwpm(experiment.dem);
            double unit = edge_weight_units(probe.flooder.graph);

            BallConfig config;
            config.T = 2.0 * unit;
            config.ball.T_max = 2.0 * unit;
            config.ball.R = 4.0 * unit;
            config.mode = options.mode;
            BallDecoder decoder = BallDecoder::from_mwpm(build_mwpm(experiment.dem), config);

            pm::Mwpm exact = build_mwpm(experiment.dem);
            pm::Mwpm cleanup = build_mwpm(experiment.dem);

            uint64_t errors_exact = 0;
            uint64_t errors_m1 = 0;
            uint64_t errors_m2 = 0;
            for (size_t s = 0; s < experiment.shots.size(); s++) {
                const auto& shot = experiment.shots[s];
                uint8_t truth = experiment.observable_flips[s];

                pm::MatchingResult stock = pm::decode_detection_events_for_up_to_64_observables(exact, shot, false);
                if ((uint8_t)(stock.obs_mask & 1) != truth)
                    errors_exact++;

                // Both front ends feed the *same* exact cleanup, so any LER difference is the front
                // end's and nothing else's.
                auto finish = [&](const HarvestResult& harvest) {
                    pm::obs_int obs = harvest.committed.obs_mask;
                    if (!harvest.residual.empty()) {
                        pm::MatchingResult rest =
                            pm::decode_detection_events_for_up_to_64_observables(cleanup, harvest.residual, false);
                        obs ^= rest.obs_mask;
                    }
                    return (uint8_t)(obs & 1);
                };
                if (finish(decoder.reference_phase1_on_g(shot)) != truth)
                    errors_m1++;
                if (finish(decoder.decode_phase1(shot)) != truth)
                    errors_m2++;
            }

            double denominator = (double)experiment.shots.size();
            double ler_exact = (double)errors_exact / denominator;
            double ler_m1 = (double)errors_m1 / denominator;
            double ler_m2 = (double)errors_m2 / denominator;
            std::printf("%4zu %8.4f %12.3e %12.3e %12.3e\n", distance, noise, ler_exact, ler_m1, ler_m2);
            emit("ler", distance, noise, 2.0, "ler_exact", ler_exact);
            emit("ler", distance, noise, 2.0, "ler_m1", ler_m1);
            emit("ler", distance, noise, 2.0, "ler_m2", ler_m2);
            emit("ler", distance, noise, 2.0, "shots", denominator);

            distances_by_p[noise].push_back((double)distance);
            ler_exact_by_p[noise].push_back(ler_exact);
            ler_m1_by_p[noise].push_back(ler_m1);
            ler_m2_by_p[noise].push_back(ler_m2);
        }
    }

    std::printf("\n=== level 3: slope of ln(LER) vs d (proportional to d_eff/d) ===\n");
    std::printf("%8s %20s %20s %20s\n", "p", "exact", "m1", "m2");
    for (const auto& [noise, distances] : distances_by_p) {
        Fit exact_fit = fit_log_ler(distances, ler_exact_by_p[noise]);
        Fit m1_fit = fit_log_ler(distances, ler_m1_by_p[noise]);
        Fit m2_fit = fit_log_ler(distances, ler_m2_by_p[noise]);
        std::printf(
            "%8.4f %10.4f+-%-8.4f %10.4f+-%-8.4f %10.4f+-%-8.4f\n",
            noise,
            exact_fit.slope,
            exact_fit.stderr_slope,
            m1_fit.slope,
            m1_fit.stderr_slope,
            m2_fit.slope,
            m2_fit.stderr_slope);
        emit("ler_fit", 0, noise, 2.0, "slope_exact", exact_fit.slope);
        emit("ler_fit", 0, noise, 2.0, "stderr_exact", exact_fit.stderr_slope);
        emit("ler_fit", 0, noise, 2.0, "slope_m1", m1_fit.slope);
        emit("ler_fit", 0, noise, 2.0, "stderr_m1", m1_fit.stderr_slope);
        emit("ler_fit", 0, noise, 2.0, "slope_m2", m2_fit.slope);
        emit("ler_fit", 0, noise, 2.0, "stderr_m2", m2_fit.stderr_slope);
    }

    if (csv.enabled)
        std::printf("\nwrote %s\n", options.csv_path.c_str());
    return 0;
}
