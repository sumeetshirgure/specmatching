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

/// The M6 exit artifact — §M6.4's benchmark list.
///
/// Five sections. The escalation and streaming budgets live in the M3 artifact, which this one does
/// not duplicate; run both.
///
///  - `latency`   — the end-to-end distribution against stock exact, per `(d, p, T)`, out to p99.99
///                  and the max, with `contaminated_shot_rate` beside every percentile, plus the
///                  `summarize()` reconciliation of the amortised mean against the measured one.
///  - `frontend`  — three curves on one axis per `(d, p)`: stock exact, M1's Phase 1 on `G`, and
///                  M2's Phase 1 on `H`. Swept down to low `p`, where the comparison is supposed to
///                  favour `H` most.
///  - `ler`       — the effective distance. Fits `LER ~ c1 * (p/p_th)^(c2 * d)` and reports
///                  `c2 = d_eff/d` for the two-phase decoder against exact MWPM. It is **expected
///                  to be equal**, since the output is exact on every shot; a departure is a Phase-1
///                  commit-policy bug, not an accuracy tradeoff. Run at `p ∈ {3e-3, 5e-3}` as well
///                  as the operating grid: at `p <= 1e-3` every LER in the campaign is zero at 10⁴
///                  shots and the fit has no signal, which has produced a vacuous Level 3 result
///                  three campaigns running.
///  - `ball`      — table cost per `d` and the §M2 structural counters, as the local-memory and
///                  interconnect budget. `BITSET` reads roughly 4x fewer bytes per defect than
///                  `SCAN`, so the row says which mode it is.
///  - `heavy`     — the top 0.1% of shots by defect count, plus every escalating shot, with the
///                  full profile breakdown.
///
/// Usage:
///   two_phase_m6_artifact [--distances 5,7,9,11,13] [--error-rates 0.0005,0.001,0.003,0.005]
///                         [--horizons 2.0] [--shots 20000] [--ler-shots 20000]
///                         [--mode scan|bitset] [--seed N] [--csv path]

#include <cstdio>
#include <map>

#include "benchmarks/two_phase/artifact_util.h"
#include "pyrematching/two_phase/driver/two_phase_decoding.h"
#include "pyrematching/two_phase/manifold/ball_decoding.h"

using namespace pm::two_phase;
using namespace pm::two_phase::artifact;

namespace {

struct Options {
    std::vector<size_t> distances = {5, 7, 9, 11, 13};
    std::vector<double> error_rates = {0.0005, 0.001, 0.003, 0.005};
    std::vector<double> horizons = {2.0};
    size_t shots = 20000;
    size_t ler_shots = 20000;
    BallGraphBuildMode mode = BallGraphBuildMode::SCAN;
    uint64_t seed = 20260804;
    std::string csv_path;
};

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

TwoPhaseConfig config_for(double horizon_multiple, double unit, BallGraphBuildMode mode) {
    TwoPhaseConfig config;
    config.T = horizon_multiple * unit;
    config.ball.T_max = horizon_multiple * unit;
    config.ball.R = 2.0 * config.ball.T_max;
    config.mode = mode;
    return config;
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

    const char* mode_name = options.mode == BallGraphBuildMode::BITSET ? "bitset" : "scan";

    // ------------------------------------------------------------- sections: latency, frontend, heavy
    std::printf(
        "\n=== end-to-end latency vs stock exact (%zu shots per point, %s build) ===\n", options.shots, mode_name);
    std::printf(
        "%4s %8s %5s %10s %10s %10s %10s %10s %9s %9s %10s\n",
        "d",
        "p",
        "T",
        "mean",
        "p50",
        "p99",
        "p999",
        "p9999",
        "stock",
        "speedup",
        "contam");
    std::vector<std::string> frontend_rows;
    std::vector<std::string> heavy_rows;

    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, distance, noise, options.shots, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double horizon_multiple : options.horizons) {
                TwoPhaseConfig config = config_for(horizon_multiple, unit, options.mode);
                config.measure_exact_reference = true;
                config.detect_preemption = true;
                auto decoder = TwoPhaseDecoder::from_detector_error_model(experiment.dem, config, NUM_DISTINCT_WEIGHTS);
                decoder.stats.keep_per_shot = true;

                std::vector<uint8_t> warmup(decoder.num_observables, 0);
                pm::total_weight_int warmup_weight = 0;
                if (!experiment.shots.empty())
                    decoder.decode_to_obs(experiment.shots[0], warmup.data(), warmup_weight, nullptr);

                std::vector<TwoPhaseProfile> profiles;
                profiles.reserve(experiment.shots.size());
                decoder.decode_batch(experiment.shots, nullptr, nullptr, true, &profiles);
                TwoPhaseSummary summary = summarize(decoder.stats);
                double mean_stock =
                    (double)decoder.stats.sum_exact_reference_ns / (double)std::max<uint64_t>(1, decoder.stats.shots);

                std::printf(
                    "%4zu %8g %5g %10.1f %10.1f %10.1f %10.1f %10.1f %9.1f %9.3f %9.4f%%\n",
                    distance,
                    noise,
                    horizon_multiple,
                    summary.measured_mean_ns,
                    summary.p50_total_ns,
                    summary.p99_total_ns,
                    summary.p999_total_ns,
                    summary.p9999_total_ns,
                    mean_stock,
                    summary.speedup_vs_stock,
                    100.0 * summary.contaminated_shot_rate);
                emit("latency", distance, noise, horizon_multiple, "shots", (double)decoder.stats.shots);
                emit(
                    "latency",
                    distance,
                    noise,
                    horizon_multiple,
                    "shots_escalated",
                    (double)decoder.stats.shots_escalated);
                emit("latency", distance, noise, horizon_multiple, "q", summary.q);
                emit("latency", distance, noise, horizon_multiple, "mean_ns", summary.measured_mean_ns);
                emit("latency", distance, noise, horizon_multiple, "p50_ns", summary.p50_total_ns);
                emit("latency", distance, noise, horizon_multiple, "p99_ns", summary.p99_total_ns);
                emit("latency", distance, noise, horizon_multiple, "p999_ns", summary.p999_total_ns);
                emit("latency", distance, noise, horizon_multiple, "p9999_ns", summary.p9999_total_ns);
                emit("latency", distance, noise, horizon_multiple, "max_ns", summary.max_total_ns);
                emit("latency", distance, noise, horizon_multiple, "stock_mean_ns", mean_stock);
                emit("latency", distance, noise, horizon_multiple, "speedup_vs_stock", summary.speedup_vs_stock);
                emit(
                    "latency",
                    distance,
                    noise,
                    horizon_multiple,
                    "contaminated_shot_rate",
                    summary.contaminated_shot_rate);
                // The M6 exit checkpoint's reconciliation, as a number rather than a judgement call.
                emit("latency", distance, noise, horizon_multiple, "amortised_mean_ns", summary.amortised_mean_ns);
                emit("latency", distance, noise, horizon_multiple, "amortisation_gap", summary.amortisation_gap);

                // --------------------------------------------------------------- heavy-shot drill-down
                //
                // Everything unusual, in one place: the top 0.1% by defect count and every shot that
                // escalated. A mean over 20 000 shots hides both.
                std::vector<size_t> order(profiles.size());
                for (size_t i = 0; i < order.size(); i++)
                    order[i] = i;
                std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                    return profiles[a].num_defects > profiles[b].num_defects;
                });
                size_t heavy_count = std::max<size_t>(1, profiles.size() / 1000);
                std::vector<size_t> selected(order.begin(), order.begin() + std::min(heavy_count, order.size()));
                for (size_t i = 0; i < profiles.size(); i++) {
                    if (profiles[i].escalated)
                        selected.push_back(i);
                }
                std::sort(selected.begin(), selected.end());
                selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
                for (size_t i : selected) {
                    const TwoPhaseProfile& profile = profiles[i];
                    std::ostringstream row;
                    row << "heavy," << distance << "," << noise << "," << horizon_multiple << ",shot_" << i << ",0";
                    heavy_rows.push_back(row.str());
                    auto heavy = [&](const std::string& key, double value) {
                        std::ostringstream heavy_row;
                        heavy_row << "heavy," << distance << "," << noise << "," << horizon_multiple << ",shot_" << i
                                  << "_" << key << "," << value;
                        heavy_rows.push_back(heavy_row.str());
                    };
                    heavy("num_defects", profile.num_defects);
                    heavy("phase1_ns", (double)profile.phase1_ns);
                    heavy("harvest_ns", (double)profile.harvest_ns);
                    heavy("escalation_ns", (double)profile.escalation_ns);
                    heavy("stock_ns", (double)profile.stock_ns);
                    heavy("total_ns", (double)profile.total_ns);
                    heavy("exact_reference_ns", (double)profile.exact_reference_ns);
                    heavy("escalated", profile.escalated ? 1 : 0);
                    heavy("contaminated", profile.contaminated ? 1 : 0);
                }

                // ------------------------------------------------------------- front-end comparison
                //
                // Three curves, measured on the same shots: stock exact on `G`, M1's Phase 1 on `G`,
                // and M2's Phase 1 on `H`. The first is a different job from the other two — it
                // returns a full matching where Phase 1 returns a partial one — so it is a separate
                // column and never enters a Phase-1-to-Phase-1 ratio.
                BallConfig ball_config;
                ball_config.T = config.T;
                ball_config.ball = config.ball;
                ball_config.mode = options.mode;
                auto ball_decoder =
                    BallDecoder::from_detector_error_model(experiment.dem, ball_config, NUM_DISTINCT_WEIGHTS);
                pm::Mwpm g_mwpm = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
                Harvester g_harvester;
                horizon_int horizon = ball_decoder.horizon;

                long long sum_h_ns = 0;
                long long sum_g_ns = 0;
                HiResTimer timer;
                size_t sampled = std::min<size_t>(experiment.shots.size(), 4000);
                for (size_t i = 0; i < sampled; i++) {
                    const auto& shot = experiment.shots[i];
                    BallProfile ball_profile;
                    ball_decoder.decode_phase1_production(shot, &ball_profile);
                    sum_h_ns += ball_profile.total_ns;

                    timer.start();
                    process_timeline_until_horizon(g_mwpm, shot, horizon);
                    g_harvester.harvest_to_obs(g_mwpm, shot);
                    sum_g_ns += timer.elapsed_ns();
                }
                double mean_h = (double)sum_h_ns / (double)std::max<size_t>(1, sampled);
                double mean_g = (double)sum_g_ns / (double)std::max<size_t>(1, sampled);
                std::ostringstream frontend;
                frontend << "  d=" << distance << " p=" << noise << " T=" << horizon_multiple
                         << "  stock=" << mean_stock << "ns  m1_on_G=" << mean_g << "ns  m2_on_H=" << mean_h
                         << "ns  m2/m1=" << (mean_g > 0 ? mean_h / mean_g : 0);
                frontend_rows.push_back(frontend.str());
                emit("frontend", distance, noise, horizon_multiple, "stock_mean_ns", mean_stock);
                emit("frontend", distance, noise, horizon_multiple, "m1_phase1_on_g_ns", mean_g);
                emit("frontend", distance, noise, horizon_multiple, "m2_phase1_on_h_ns", mean_h);
                emit("frontend", distance, noise, horizon_multiple, "speedup_vs_m1", mean_h > 0 ? mean_g / mean_h : 0);
                emit("frontend", distance, noise, horizon_multiple, "shots_sampled", (double)sampled);
            }
        }
    }

    std::printf("\n=== front-end comparison (stock exact vs M1 Phase 1 on G vs M2 Phase 1 on H) ===\n");
    for (const std::string& row : frontend_rows)
        std::printf("%s\n", row.c_str());

    // ---------------------------------------------------------------------------- section: ball
    std::printf("\n=== ball table cost and structural counters (%s build) ===\n", mode_name);
    std::printf(
        "%4s %10s %12s %12s %10s %10s %12s %12s %12s\n",
        "d",
        "nodes",
        "entries",
        "bytes",
        "mean|B|",
        "max|B|",
        "isect_B/def",
        "hbld_edges",
        "mwpm_init");
    for (size_t distance : options.distances) {
        double noise = options.error_rates.back();
        Experiment experiment =
            generate(distance, distance, noise, std::min<size_t>(options.shots, 2000), options.seed);
        pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
        double unit = edge_weight_units(probe.flooder.graph);

        BallConfig ball_config;
        ball_config.T = options.horizons.front() * unit;
        ball_config.ball.T_max = ball_config.T;
        ball_config.ball.R = 2.0 * ball_config.ball.T_max;
        ball_config.mode = options.mode;
        // The counters are structural — a function of the shot and the tables, not of when they
        // were measured — so collecting them in this separate, untimed pass costs nothing in
        // fidelity and keeps the byte counting out of the timing table above.
        ball_config.collect_structural_counters = true;
        auto ball_decoder = BallDecoder::from_detector_error_model(experiment.dem, ball_config, NUM_DISTINCT_WEIGHTS);

        BallAggregateStats aggregate;
        for (const auto& shot : experiment.shots) {
            BallProfile profile;
            Phase1Outcome outcome = ball_decoder.decode_phase1_production(shot, &profile);
            aggregate.accumulate(profile, outcome.harvest);
        }
        BallSummary ball_summary = summarize_ball(aggregate);
        const BallStats& table_stats = ball_decoder.tables.stats;

        std::printf(
            "%4zu %10llu %12llu %12llu %10.1f %10llu %12.1f %12.1f %12.1f\n",
            distance,
            (unsigned long long)table_stats.num_nodes,
            (unsigned long long)table_stats.total_entries,
            (unsigned long long)table_stats.bytes_total,
            table_stats.mean_ball_size,
            (unsigned long long)table_stats.max_ball_size,
            ball_summary.isect_bytes_per_defect,
            ball_summary.mean_hbld_edges_written,
            ball_summary.mean_mwpm_init_elements);
        emit("ball", distance, noise, 0, "bytes_total", (double)table_stats.bytes_total);
        emit("ball", distance, noise, 0, "bytes_paths", (double)table_stats.bytes_paths);
        emit("ball", distance, noise, 0, "mean_ball_size", table_stats.mean_ball_size);
        emit("ball", distance, noise, 0, "max_ball_size", (double)table_stats.max_ball_size);
        emit("ball", distance, noise, 0, "mean_degree_at_2T", ball_summary.mean_degree);
        emit("ball", distance, noise, 0, "max_degree_at_2T", ball_summary.max_degree);
        emit("ball", distance, noise, 0, "mode_bitset", ball_summary.mode_bitset ? 1 : 0);
        emit("ball", distance, noise, 0, "isect_bytes_per_defect", ball_summary.isect_bytes_per_defect);
        emit(
            "ball",
            distance,
            noise,
            0,
            "isect_bytes_other_mode_per_defect",
            ball_summary.isect_bytes_other_mode_per_defect);
        emit("ball", distance, noise, 0, "hbld_edges_written", ball_summary.mean_hbld_edges_written);
        emit("ball", distance, noise, 0, "mwpm_init_elements", ball_summary.mean_mwpm_init_elements);
        emit("ball", distance, noise, 0, "compile_wall_seconds", table_stats.compile_wall_seconds);
    }

    // ----------------------------------------------------------------------------- section: ler
    //
    // The effective distance. Both decoders are exact MWPM, so `c2 = d_eff/d` is expected to be
    // *equal*, and the point of measuring it is that a departure is a commit-policy bug. Run at the
    // high-`p` end as well: at `p <= 1e-3` every LER is zero at these shot counts and the fit has
    // no signal at all.
    std::printf("\n=== effective distance (%zu shots per point) ===\n", options.ler_shots);
    std::printf(
        "%8s %10s %16s %16s %8s %10s\n", "p", "d_points", "slope_two_phase", "slope_exact", "agree", "disagree");
    std::vector<double> ler_error_rates = options.error_rates;
    for (double extra : {0.003, 0.005}) {
        if (std::find(ler_error_rates.begin(), ler_error_rates.end(), extra) == ler_error_rates.end())
            ler_error_rates.push_back(extra);
    }
    for (double noise : ler_error_rates) {
        std::vector<double> distances;
        std::vector<double> two_phase_lers;
        std::vector<double> exact_lers;
        // The two decoders make the *same prediction on every shot* — that is what M3 buys, and it
        // is asserted bit-exactly by §M4.2's test. Counting the disagreements here says so directly;
        // the fit below can only ever be a weaker, statistical restatement of this number, so it is
        // the one to read first.
        size_t prediction_disagreements = 0;
        for (size_t distance : options.distances) {
            Experiment experiment = generate(distance, distance, noise, options.ler_shots, options.seed + 1);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);
            TwoPhaseConfig config = config_for(options.horizons.front(), unit, options.mode);
            auto decoder = TwoPhaseDecoder::from_detector_error_model(experiment.dem, config, NUM_DISTINCT_WEIGHTS);
            pm::Mwpm exact_mwpm = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);

            size_t two_phase_errors = 0;
            size_t exact_errors = 0;
            std::vector<uint8_t> obs(decoder.num_observables, 0);
            for (size_t i = 0; i < experiment.shots.size(); i++) {
                pm::total_weight_int weight = 0;
                decoder.decode_to_obs(experiment.shots[i], obs.data(), weight, nullptr);
                if ((obs[0] & 1) != experiment.observable_flips[i])
                    two_phase_errors++;

                pm::MatchingResult exact =
                    pm::decode_detection_events_for_up_to_64_observables(exact_mwpm, experiment.shots[i], false);
                if ((exact.obs_mask & 1) != experiment.observable_flips[i])
                    exact_errors++;
                if ((obs[0] & 1) != (exact.obs_mask & 1))
                    prediction_disagreements++;
            }
            double shots = (double)experiment.shots.size();
            distances.push_back((double)distance);
            two_phase_lers.push_back((double)two_phase_errors / shots);
            exact_lers.push_back((double)exact_errors / shots);
            emit("ler", distance, noise, options.horizons.front(), "two_phase_ler", (double)two_phase_errors / shots);
            emit("ler", distance, noise, options.horizons.front(), "exact_ler", (double)exact_errors / shots);
            emit("ler", distance, noise, options.horizons.front(), "shots", shots);
        }
        Fit two_phase_fit = fit_log_ler(distances, two_phase_lers);
        Fit exact_fit = fit_log_ler(distances, exact_lers);
        bool agree = two_phase_fit.points >= 3 && exact_fit.points >= 3 &&
                     std::abs(two_phase_fit.slope - exact_fit.slope) <=
                         2.0 * std::sqrt(
                                   two_phase_fit.stderr_slope * two_phase_fit.stderr_slope +
                                   exact_fit.stderr_slope * exact_fit.stderr_slope);
        std::printf(
            "%8g %10zu %10.4f+-%.4f %10.4f+-%.4f %8s %10zu\n",
            noise,
            two_phase_fit.points,
            two_phase_fit.slope,
            two_phase_fit.stderr_slope,
            exact_fit.slope,
            exact_fit.stderr_slope,
            two_phase_fit.points < 3 ? "no signal" : (agree ? "yes" : "NO"),
            prediction_disagreements);
        // The slope of `ln(LER)` against `d` at fixed `p`. `c2 = d_eff/d` is *proportional* to it —
        // the constant is `ln(p / p_th)` — so the two front ends are compared on the slope, which is
        // what has to agree, rather than on an absolute `c2` that would need a threshold estimate.
        emit("ler", 0, noise, options.horizons.front(), "slope_two_phase", two_phase_fit.slope);
        emit("ler", 0, noise, options.horizons.front(), "slope_two_phase_stderr", two_phase_fit.stderr_slope);
        emit("ler", 0, noise, options.horizons.front(), "slope_exact", exact_fit.slope);
        emit("ler", 0, noise, options.horizons.front(), "slope_exact_stderr", exact_fit.stderr_slope);
        emit("ler", 0, noise, options.horizons.front(), "fit_points", (double)two_phase_fit.points);
        emit("ler", 0, noise, options.horizons.front(), "prediction_disagreements", (double)prediction_disagreements);
    }

    for (const std::string& row : heavy_rows)
        csv.row(row);
    if (csv.enabled)
        std::printf("\nwrote %s\n", options.csv_path.c_str());
    return 0;
}
