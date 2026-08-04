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

/// The M3 exit artifact — the escalation budget.
///
/// Four sections, matching the four measured lines of the M3 exit checkpoint:
///
///  - `q`         — the escalation rate per `(d, p, T)`, **always with its shot count**. A `q` of
///                  zero below `1 / shots` is the resolution floor of the campaign, not a
///                  measurement: §M3.0's own campaign reported `q ∈ {0, 0.0005}` at 2000 shots and
///                  both were artifacts of `1/2000`. Size `--shots` accordingly; 10⁶ at the
///                  operating point, not 2000.
///  - `cost`      — `C_phase1`, `C_escalation`, and `mean_stock_ns_on_escalated`, which is §M3.2's
///                  `stock` term measured on precisely the shots that escalate rather than
///                  estimated as a mean over all of them. The escalation path produces it for free.
///  - `latency`   — the end-to-end distribution with escalation live, out to **p99.99 and the max**.
///                  At `q ~ 3e-4` the spike sits near p99.97, so p99 and p999 do not show it and
///                  are not sufficient evidence. `contaminated_shot_rate` sits beside every
///                  percentile; without it a percentile is a measurement of the scheduler.
///  - `streaming` — worst-case escalation latency, and the input buffer depth that absorbs it at
///                  the target round rate. This is the number a real-time deployment needs and it
///                  does not follow from any mean.
///
/// Horizons are quoted in multiples of one lattice edge weight.
///
/// Usage:
///   two_phase_m3_artifact [--distances 11,13] [--error-rates 0.0005,0.001]
///                         [--horizons 1.5,2.0] [--shots 200000] [--chunk 50000] [--rounds 0]
///                         [--round-ns 1000] [--mode scan|bitset] [--seed N] [--csv path]
///
/// `--rounds 0` means "as many rounds as the distance", which is the `d`-round shot the M1 campaign
/// used. `q` is extensive in spacetime volume (M1 result 3), so a deployment's `q` must be sized at
/// its intended `(d, rounds)` and not read off this table.

#include <cstdio>

#include "benchmarks/two_phase/artifact_util.h"
#include "pyrematching/two_phase/driver/two_phase_decoding.h"

using namespace pm::two_phase;
using namespace pm::two_phase::artifact;

namespace {

struct Options {
    std::vector<size_t> distances = {11, 13};
    std::vector<double> error_rates = {0.0005, 0.001};
    std::vector<double> horizons = {1.5, 2.0};
    size_t shots = 200000;
    /// Shots sampled and decoded at a time. Bounds the harness's memory rather than the campaign's
    /// length, so a 10⁷-shot budget run costs time and not a gigabyte of resident syndromes.
    size_t chunk = 50000;
    /// 0 means "rounds = distance".
    size_t rounds = 0;
    /// The target syndrome-extraction round time, in nanoseconds. The streaming budget is quoted
    /// against it; 1 µs is a plausible superconducting round.
    double round_ns = 1000;
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
        } else if (flag == "--chunk") {
            options.chunk = std::max<size_t>(1, std::stoul(next()));
        } else if (flag == "--rounds") {
            options.rounds = std::stoul(next());
        } else if (flag == "--round-ns") {
            options.round_ns = std::stod(next());
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

struct Point {
    size_t distance{0};
    size_t rounds{0};
    double error_rate{0};
    double horizon_multiple{0};
    TwoPhaseAggregateStats stats;
    TwoPhaseSummary summary;
    /// Stock exact decode on the same shots, as a paired measurement.
    double mean_stock_ns{0};
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
    csv.row("section,d,rounds,p,T,key,value");
    auto emit = [&](const std::string& section, const Point& point, const std::string& key, double value) {
        std::ostringstream row;
        row << section << "," << point.distance << "," << point.rounds << "," << point.error_rate << ","
            << point.horizon_multiple << "," << key << "," << value;
        csv.row(row.str());
    };

    std::vector<Point> points;
    for (size_t distance : options.distances) {
        size_t rounds = options.rounds == 0 ? distance : options.rounds;
        for (double noise : options.error_rates) {
            ShotSampler sampler = ShotSampler::make(distance, rounds, noise, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double horizon_multiple : options.horizons) {
                Point point;
                point.distance = distance;
                point.rounds = rounds;
                point.error_rate = noise;
                point.horizon_multiple = horizon_multiple;

                TwoPhaseConfig config;
                config.T = horizon_multiple * unit;
                config.ball.T_max = horizon_multiple * unit;
                config.ball.R = 2.0 * config.ball.T_max;
                config.mode = options.mode;
                // Paired against stock on the same shot, and marked contaminated when the
                // scheduler intervened (§M6.4 timer discipline).
                config.measure_exact_reference = true;
                config.detect_preemption = true;

                auto decoder = TwoPhaseDecoder::from_detector_error_model(sampler.dem, config, NUM_DISTINCT_WEIGHTS);
                decoder.stats.keep_per_shot = true;

                // Shots are sampled and decoded in chunks: the campaigns this section is for are
                // 10⁶ and 10⁷ shots, and materialising those costs a gigabyte before the decoder
                // has allocated anything. The stream depends on the seed alone, not on the chunk
                // size. Each horizon re-seeds, so every point in the sweep sees the same shots.
                ShotSampler point_sampler = ShotSampler::make(distance, rounds, noise, options.seed);
                bool warmed_up = false;
                std::vector<uint8_t> obs_scratch;
                for (size_t done = 0; done < options.shots;) {
                    size_t take = std::min(options.chunk, options.shots - done);
                    point_sampler.sample(take);
                    if (!warmed_up && !point_sampler.shots.empty()) {
                        // One warmup shot, so that the arena's steady state — not its first
                        // allocation — is what the percentiles see.
                        obs_scratch.assign(decoder.num_observables, 0);
                        pm::total_weight_int weight = 0;
                        decoder.decode_to_obs(point_sampler.shots[0], obs_scratch.data(), weight, nullptr);
                        warmed_up = true;
                    }
                    decoder.decode_batch(point_sampler.shots, nullptr, nullptr, true, nullptr);
                    done += take;
                }

                point.stats = decoder.stats;
                point.summary = summarize(point.stats);
                point.mean_stock_ns =
                    (double)point.stats.sum_exact_reference_ns / (double)std::max<uint64_t>(1, point.stats.shots);
                points.push_back(std::move(point));
                std::fprintf(
                    stderr,
                    "  d=%zu rounds=%zu p=%g T=%g done (%llu shots, %llu escalated)\n",
                    distance,
                    rounds,
                    noise,
                    horizon_multiple,
                    (unsigned long long)points.back().stats.shots,
                    (unsigned long long)points.back().stats.shots_escalated);
            }
        }
    }

    // ------------------------------------------------------------------------------- section: q
    std::printf("\n=== escalation rate (always with its shot count) ===\n");
    std::printf(
        "%4s %7s %8s %5s %10s %10s %12s %14s\n", "d", "rounds", "p", "T", "shots", "escalated", "q", "resolution");
    for (const Point& point : points) {
        double resolution = 1.0 / (double)std::max<uint64_t>(1, point.stats.shots);
        std::printf(
            "%4zu %7zu %8g %5g %10llu %10llu %12.3e %14.3e%s\n",
            point.distance,
            point.rounds,
            point.error_rate,
            point.horizon_multiple,
            (unsigned long long)point.stats.shots,
            (unsigned long long)point.stats.shots_escalated,
            point.summary.q,
            resolution,
            point.stats.shots_escalated == 0 ? "  <- at or below the floor" : "");
        emit("q", point, "shots", (double)point.stats.shots);
        emit("q", point, "shots_escalated", (double)point.stats.shots_escalated);
        emit("q", point, "q", point.summary.q);
        emit("q", point, "resolution_floor", resolution);
        emit("q", point, "mean_residual_density", point.summary.mean_residual_density);
    }

    // ---------------------------------------------------------------------------- section: cost
    std::printf("\n=== escalation cost (ns) ===\n");
    std::printf(
        "%4s %8s %5s %11s %13s %13s %11s %11s %11s\n",
        "d",
        "p",
        "T",
        "C_phase1",
        "C_escalation",
        "stock_on_esc",
        "mean_stock",
        "amortised",
        "penalty");
    for (const Point& point : points) {
        std::printf(
            "%4zu %8g %5g %11.1f %13.1f %13.1f %11.1f %11.1f %10.4f%%\n",
            point.distance,
            point.error_rate,
            point.horizon_multiple,
            point.summary.c_phase1,
            point.summary.c_escalation,
            point.summary.mean_stock_ns_on_escalated,
            point.mean_stock_ns,
            point.summary.amortised_mean_ns,
            100.0 * point.summary.amortised_penalty);
        emit("cost", point, "c_phase1_ns", point.summary.c_phase1);
        emit("cost", point, "c_escalation_ns", point.summary.c_escalation);
        // §M3.2 asked for this to replace its estimate. `mean_stock_ns` is that estimate — the mean
        // over *all* shots — and it is printed next to it so the gap is visible rather than
        // asserted: escalating shots have surviving trees at `T` and are harder than average.
        emit("cost", point, "mean_stock_ns_on_escalated", point.summary.mean_stock_ns_on_escalated);
        emit("cost", point, "mean_stock_ns_all_shots", point.mean_stock_ns);
        emit("cost", point, "amortised_mean_ns", point.summary.amortised_mean_ns);
        emit("cost", point, "measured_mean_ns", point.summary.measured_mean_ns);
        emit("cost", point, "amortisation_gap", point.summary.amortisation_gap);
        emit("cost", point, "amortised_penalty", point.summary.amortised_penalty);
        emit("cost", point, "escalation_cost_ratio", point.summary.escalation_cost_ratio);
    }

    // ------------------------------------------------------------------------- section: latency
    std::printf("\n=== end-to-end latency with escalation live (ns) ===\n");
    std::printf(
        "%4s %8s %5s %10s %10s %10s %10s %10s %10s %12s\n",
        "d",
        "p",
        "T",
        "mean",
        "p50",
        "p99",
        "p999",
        "p9999",
        "max",
        "contam");
    for (const Point& point : points) {
        std::printf(
            "%4zu %8g %5g %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %11.4f%%\n",
            point.distance,
            point.error_rate,
            point.horizon_multiple,
            point.summary.measured_mean_ns,
            point.summary.p50_total_ns,
            point.summary.p99_total_ns,
            point.summary.p999_total_ns,
            point.summary.p9999_total_ns,
            point.summary.max_total_ns,
            100.0 * point.summary.contaminated_shot_rate);
        emit("latency", point, "mean_ns", point.summary.measured_mean_ns);
        emit("latency", point, "p50_ns", point.summary.p50_total_ns);
        emit("latency", point, "p99_ns", point.summary.p99_total_ns);
        emit("latency", point, "p999_ns", point.summary.p999_total_ns);
        emit("latency", point, "p9999_ns", point.summary.p9999_total_ns);
        emit("latency", point, "max_ns", point.summary.max_total_ns);
        emit("latency", point, "contaminated_shot_rate", point.summary.contaminated_shot_rate);
        emit("latency", point, "speedup_vs_stock", point.summary.speedup_vs_stock);
        emit("latency", point, "stock_mean_ns", point.mean_stock_ns);
        emit("latency", point, "stock_p99_ns", percentile(point.stats.per_shot_exact_ns, 0.99));
        emit("latency", point, "stock_p9999_ns", percentile(point.stats.per_shot_exact_ns, 0.9999));
    }

    // ----------------------------------------------------------------------- section: streaming
    //
    // The buffer has to absorb the worst escalating shot while the source keeps producing rounds.
    // With `rounds` rounds arriving every `round_ns` each, a shot's budget is `rounds * round_ns`;
    // anything an escalating shot spends beyond that has to be soaked up by queued input, and the
    // depth in shots is the overrun divided by the per-shot budget.
    std::printf("\n=== streaming budget (round time %.0f ns) ===\n", options.round_ns);
    std::printf(
        "%4s %8s %5s %13s %13s %13s %12s %10s\n",
        "d",
        "p",
        "T",
        "shot_budget",
        "worst_esc",
        "worst_any",
        "overrun",
        "buffer");
    for (const Point& point : points) {
        double shot_budget_ns = (double)point.rounds * options.round_ns;
        double worst_escalated = point.summary.max_escalated_total_ns;
        double worst_any = point.summary.max_total_ns;
        double overrun = std::max(0.0, worst_any - shot_budget_ns);
        double buffer_depth = shot_budget_ns > 0 ? std::ceil(overrun / shot_budget_ns) : 0;
        std::printf(
            "%4zu %8g %5g %13.1f %13.1f %13.1f %12.1f %10.0f\n",
            point.distance,
            point.error_rate,
            point.horizon_multiple,
            shot_budget_ns,
            worst_escalated,
            worst_any,
            overrun,
            buffer_depth);
        emit("streaming", point, "round_ns", options.round_ns);
        emit("streaming", point, "shot_budget_ns", shot_budget_ns);
        emit("streaming", point, "worst_escalated_ns", worst_escalated);
        emit("streaming", point, "worst_any_ns", worst_any);
        emit("streaming", point, "overrun_ns", overrun);
        emit("streaming", point, "buffer_depth_shots", buffer_depth);
    }

    if (csv.enabled)
        std::printf("\nwrote %s\n", options.csv_path.c_str());
    return 0;
}
