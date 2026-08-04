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

/// The M2.9 exit artifact.
///
/// Four sections, matching the M2.9 exit checkpoint. Unlike the M1 and M2 artifacts, the headline
/// numbers here are **counts and depths, not nanoseconds**: the governing metric from the M2 exit
/// read onwards is the critical path in dependent operations, and a nanosecond figure is a
/// statement about this machine's cache hierarchy rather than about the architecture the design is
/// aimed at. Wall time appears only where the design explicitly asks for it — the software A/B of
/// §M2.9.6, which asks how much of harvest's share of the total was enumeration.
///
///  - `measure` — §M2.9.6 measurements 1–4 per `(d, p, T)`: the `largest_tree_size` distribution,
///                blossom nesting depth and member count (split into all blossoms and the exposed
///                root blossoms that §M2.9.4's base descent actually walks), blossom formations per
///                shot, and the dependent-event chain depth of the solve. Plus harvest's own
///                modelled depth, and the ratio the exit checkpoint asks for.
///  - `ab`      — the software A/B: `harvest_ns` split into enumerate / base descent / shatter /
///                reduce, and `speedup_vs_m1` re-measured with the restructured harvest so that the
///                M2 exit number is comparable.
///  - `split`   — §M2.9.5 E1 and E2, over the full corpus, per horizon. **These are the gate**: the
///                deferred shatter is not to be implemented unless both pass.
///  - `depth`   — the exposed-root-blossom nesting depth histogram, which is the measurement the
///                exit checkpoint requires be written down next to §M2.9.4's decision.
///
/// Horizons are quoted in multiples of one lattice edge weight, the normalisation under which the
/// M1 results are `d`-independent.
///
/// Usage:
///   two_phase_m29_artifact [--distances 5,7,9,11,13] [--error-rates 0.001,0.003,0.005]
///                          [--horizons 1.0,1.5,2.0] [--shots 2000] [--seed N] [--csv path]

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "pyrematching/sparse_blossom/driver/mwpm_decoding.h"
#include "pyrematching/two_phase/manifold/ball_decoding.h"
#include "pyrematching/two_phase/perf/ball_profile.h"
#include "pyrematching/two_phase/truncation/split_experiment.h"
#include "stim.h"

using namespace pm::two_phase;

namespace {

const pm::weight_int NUM_DISTINCT_WEIGHTS = 10001;

struct Options {
    std::vector<size_t> distances = {5, 7, 9, 11, 13};
    std::vector<double> error_rates = {0.001, 0.003, 0.005};
    std::vector<double> horizons = {1.0, 1.5, 2.0};
    size_t shots = 2000;
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

    Experiment experiment;
    experiment.dem = stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false);
    for (size_t shot = 0; shot < num_shots; shot++) {
        experiment.shots.emplace_back();
        for (size_t d = 0; d < num_detectors; d++) {
            if (dets[d][shot])
                experiment.shots.back().push_back(d);
        }
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

    // ------------------------------------------------------------- section: measurements 1–4 + A/B
    std::printf("\n=== M2.9.6 measurements, in counts and depths (%zu shots per point) ===\n", options.shots);
    std::printf(
        "%4s %8s %5s %8s %8s %8s %8s %8s %8s %8s %8s %8s %7s\n",
        "d",
        "p",
        "T",
        "tree_max",
        "form/sh",
        "nest_max",
        "memb_max",
        "xnest",
        "xmemb",
        "solve_d",
        "harv_d",
        "events",
        "share");
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, noise, options.shots, options.seed + distance * 131);
            pm::Mwpm probe = build_mwpm(experiment.dem);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double t_edges : options.horizons) {
                BallConfig config;
                config.T = t_edges * unit;
                config.ball.T_max = std::max(2.0, t_edges) * unit;
                config.ball.R = 2.0 * config.ball.T_max;
                config.collect_harvest_diagnostics = true;
                BallDecoder decoder = BallDecoder::from_mwpm(build_mwpm(experiment.dem), config);

                BallAggregateStats stats;
                BallProfile profile;
                if (!experiment.shots.empty())
                    decoder.decode_phase1(experiment.shots[0], &profile);
                for (const auto& shot : experiment.shots) {
                    HarvestResult harvest = decoder.decode_phase1(shot, &profile);
                    stats.accumulate(profile, harvest);
                }
                BallSummary summary = summarize_ball(stats);

                std::printf(
                    "%4zu %8.4f %5.2f %8.0f %8.2f %8.0f %8.0f %8.0f %8.0f %8.1f %8.1f %8.1f %7.3f\n",
                    distance,
                    noise,
                    t_edges,
                    summary.max_largest_tree_size,
                    summary.mean_blossom_formations,
                    summary.max_blossom_nesting_depth,
                    (double)stats.max_blossom_members,
                    (double)stats.max_exposed_blossom_depth,
                    (double)stats.max_exposed_blossom_members,
                    summary.mean_solve_dependent_depth,
                    summary.mean_harvest_dependent_depth,
                    summary.mean_solve_events,
                    summary.harvest_share_of_critical_path);

                emit("measure", distance, noise, t_edges, "max_largest_tree_size", summary.max_largest_tree_size);
                emit("measure", distance, noise, t_edges, "mean_largest_tree_size", summary.mean_largest_tree_size);
                emit("measure", distance, noise, t_edges, "mean_blossom_formations", summary.mean_blossom_formations);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "mean_matched_blossom_shatters",
                    summary.mean_matched_blossom_shatters);
                emit("measure", distance, noise, t_edges, "max_blossom_nesting_depth", summary.max_blossom_nesting_depth);
                emit("measure", distance, noise, t_edges, "max_blossom_members", (double)stats.max_blossom_members);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "max_exposed_blossom_depth",
                    (double)stats.max_exposed_blossom_depth);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "max_exposed_blossom_members",
                    (double)stats.max_exposed_blossom_members);
                emit(
                    "measure", distance, noise, t_edges, "mean_solve_dependent_depth", summary.mean_solve_dependent_depth);
                emit(
                    "measure", distance, noise, t_edges, "max_solve_dependent_depth", summary.max_solve_dependent_depth);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "mean_harvest_dependent_depth",
                    summary.mean_harvest_dependent_depth);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "max_harvest_dependent_depth",
                    summary.max_harvest_dependent_depth);
                emit("measure", distance, noise, t_edges, "mean_solve_events", summary.mean_solve_events);
                emit(
                    "measure",
                    distance,
                    noise,
                    t_edges,
                    "harvest_share_of_critical_path",
                    summary.harvest_share_of_critical_path);
                for (size_t bin = 0; bin < BallAggregateStats::DEPTH_HIST_BINS; bin++) {
                    emit(
                        "depth",
                        distance,
                        noise,
                        t_edges,
                        "exposed_depth_" + std::to_string(bin),
                        (double)stats.exposed_depth_hist[bin]);
                }

            }
        }
    }

    // -------------------------------------------------------------------------- section: the A/B
    //
    // The software A/B §M2.9.6 asks for, run as two passes of one process so that the difference —
    // a few percent of a shot — is not swamped by run-to-run noise. Both passes decode the same
    // shots on the same graph with the same horizon and differ in one flag, and H1 says they
    // produce identical output; all that moves is how the committed set is enumerated.
    std::printf("\n=== M2.9.6 A/B: M1.3's enumeration vs §M2.9.1's, same process, same shots ===\n");
    std::printf(
        "%4s %8s %5s %10s %10s %8s %9s %9s %9s %9s %9s\n",
        "d",
        "p",
        "T",
        "harv_m1_ns",
        "harv_new_ns",
        "ratio",
        "frac_m1",
        "frac_new",
        "enumerate",
        "descent",
        "shatter");
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, noise, options.shots, options.seed + distance * 131);
            pm::Mwpm probe = build_mwpm(experiment.dem);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double t_edges : options.horizons) {
                BallConfig base_config;
                base_config.T = t_edges * unit;
                base_config.ball.T_max = std::max(2.0, t_edges) * unit;
                base_config.ball.R = 2.0 * base_config.ball.T_max;

                auto run = [&](bool legacy, bool diagnostics) {
                    BallConfig config = base_config;
                    config.use_legacy_harvest_enumeration = legacy;
                    config.collect_harvest_diagnostics = diagnostics;
                    BallDecoder decoder = BallDecoder::from_mwpm(build_mwpm(experiment.dem), config);
                    BallAggregateStats stats;
                    BallProfile profile;
                    // One warmup shot, so neither pass is timed with its arenas still growing.
                    if (!experiment.shots.empty())
                        decoder.decode_phase1(experiment.shots[0], &profile);
                    for (const auto& shot : experiment.shots) {
                        HarvestResult harvest = decoder.decode_phase1(shot, &profile);
                        stats.accumulate(profile, harvest);
                    }
                    return stats;
                };

                // Alternate the order across horizons so that neither flavour systematically gets
                // the colder cache.
                bool legacy_first = ((size_t)(t_edges * 100) % 2) == 0;
                BallAggregateStats legacy_stats;
                BallAggregateStats new_stats;
                if (legacy_first) {
                    legacy_stats = run(true, false);
                    new_stats = run(false, false);
                } else {
                    new_stats = run(false, false);
                    legacy_stats = run(true, false);
                }
                BallSummary legacy_summary = summarize_ball(legacy_stats);
                BallSummary new_summary = summarize_ball(new_stats);

                // A third pass, untimed against the two above, purely for the stage split.
                BallSummary split_summary = summarize_ball(run(false, true));

                double legacy_harvest_ns = (double)legacy_stats.sum_harvest_ns / (double)std::max<uint64_t>(1, legacy_stats.shots);
                double new_harvest_ns = (double)new_stats.sum_harvest_ns / (double)std::max<uint64_t>(1, new_stats.shots);
                double ratio = new_harvest_ns > 0 ? legacy_harvest_ns / new_harvest_ns : 0.0;

                std::printf(
                    "%4zu %8.4f %5.2f %10.0f %10.0f %8.2f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
                    distance,
                    noise,
                    t_edges,
                    legacy_harvest_ns,
                    new_harvest_ns,
                    ratio,
                    legacy_summary.frac_harvest,
                    new_summary.frac_harvest,
                    split_summary.frac_harvest_enumerate,
                    split_summary.frac_harvest_base_descent,
                    split_summary.frac_harvest_shatter);

                emit("ab", distance, noise, t_edges, "harvest_ns_legacy", legacy_harvest_ns);
                emit("ab", distance, noise, t_edges, "harvest_ns_direct", new_harvest_ns);
                emit("ab", distance, noise, t_edges, "harvest_speedup", ratio);
                emit("ab", distance, noise, t_edges, "frac_harvest_legacy", legacy_summary.frac_harvest);
                emit("ab", distance, noise, t_edges, "frac_harvest_direct", new_summary.frac_harvest);
                emit("ab", distance, noise, t_edges, "mean_total_ns_legacy", legacy_summary.mean_total_ns);
                emit("ab", distance, noise, t_edges, "mean_total_ns_direct", new_summary.mean_total_ns);
                emit(
                    "ab",
                    distance,
                    noise,
                    t_edges,
                    "frac_harvest_enumerate",
                    split_summary.frac_harvest_enumerate);
                emit(
                    "ab",
                    distance,
                    noise,
                    t_edges,
                    "frac_harvest_base_descent",
                    split_summary.frac_harvest_base_descent);
                emit("ab", distance, noise, t_edges, "frac_harvest_shatter", split_summary.frac_harvest_shatter);
                emit("ab", distance, noise, t_edges, "frac_harvest_reduce", split_summary.frac_harvest_reduce);
            }
        }
    }

    // ------------------------------------------------------------------ section: §M2.9.5 E1 and E2
    std::printf("\n=== M2.9.5 E1/E2: does the committed obs / weight depend on the split point? ===\n");
    std::printf(
        "%4s %8s %5s %10s %10s %10s %10s %10s %10s %10s\n",
        "d",
        "p",
        "T",
        "matches",
        "e1_moved",
        "e1_rate",
        "interior",
        "e2_moved",
        "cyc_obs!=0",
        "max_nest");
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, noise, options.shots, options.seed + distance * 131);
            pm::Mwpm probe = build_mwpm(experiment.dem);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double t_edges : options.horizons) {
                pm::Mwpm mwpm = build_mwpm(experiment.dem);
                Harvester harvester;
                auto horizon = to_time_units(t_edges * unit, mwpm.flooder.graph.normalising_constant);

                SplitExperimentStats stats;
                for (const auto& shot : experiment.shots) {
                    process_timeline_until_horizon(mwpm, shot, horizon);
                    probe_split_dependence(mwpm, stats);
                    harvester.harvest_to_obs(mwpm, shot);
                }

                double rate = stats.splittable_matches
                                  ? (double)stats.obs_differed / (double)stats.splittable_matches
                                  : 0.0;
                std::printf(
                    "%4zu %8.4f %5.2f %10llu %10llu %10.5f %10llu %10llu %6llu/%-4llu %10d\n",
                    distance,
                    noise,
                    t_edges,
                    (unsigned long long)stats.splittable_matches,
                    (unsigned long long)stats.obs_differed,
                    rate,
                    (unsigned long long)stats.interior_obs_nonzero,
                    (unsigned long long)stats.weight_differed,
                    (unsigned long long)stats.cycles_with_nonzero_obs,
                    (unsigned long long)stats.cycles_examined,
                    stats.max_nesting_depth);

                emit("split", distance, noise, t_edges, "interior_obs_nonzero", (double)stats.interior_obs_nonzero);
                emit("split", distance, noise, t_edges, "cycles_examined", (double)stats.cycles_examined);
                emit(
                    "split",
                    distance,
                    noise,
                    t_edges,
                    "cycles_with_nonzero_obs",
                    (double)stats.cycles_with_nonzero_obs);
                emit("split", distance, noise, t_edges, "splittable_matches", (double)stats.splittable_matches);
                emit(
                    "split",
                    distance,
                    noise,
                    t_edges,
                    "shots_with_matched_blossoms",
                    (double)stats.shots_with_matched_blossoms);
                emit("split", distance, noise, t_edges, "e1_obs_differed", (double)stats.obs_differed);
                emit("split", distance, noise, t_edges, "e1_rate", rate);
                emit("split", distance, noise, t_edges, "e2_weight_differed", (double)stats.weight_differed);
                emit("split", distance, noise, t_edges, "max_nesting_depth", (double)stats.max_nesting_depth);
            }
        }
    }

    if (csv.enabled)
        std::printf("\nwrote %s\n", options.csv_path.c_str());
    return 0;
}
