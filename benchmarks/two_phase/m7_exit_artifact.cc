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

/// The M7 exit artifact — §M7.8's measurements.
///
/// Two sections, in the order the exit checkpoint reads them.
///
///  - `identity`   — §M7.6 level 1 and level 3, as a **per-shot disagreement count** against stock
///                   exact decode on `G`. Expect 0. This is read instead of an LER fit: the output
///                   is exact MWPM on every shot, so a fit would only re-measure exact MWPM's own
///                   `d_eff/d`, and at `p <= 1e-3` every LER in a 10⁴-shot campaign is zero anyway
///                   (§M6.4's vacuous-fit warning). The count has signal at every `p`.
///  - `speedup`    — **the read.** Stock sparse blossom on `G` against §M2's critical path of the
///                   same stock solver on `H`, on identical shots, per `(d, p, T)`:
///
///                       speedup = stock_ns(G) / (blossom_on_h_ns + dual_scan_ns + harvest_ns)
///
///                   Only one front end runs. Both sides are *stock* blossom — untruncated, no
///                   horizon in the loop — so the ratio is the ball graph's own win and nothing
///                   else: the same solver on a smaller problem. The landed truncated-`H` scheme is
///                   not built, not timed and not compared against; §M7 replaces its front end
///                   rather than racing it.
///
///                   The denominator is §M2's critical path — the stages that sit between a syndrome
///                   arriving and a correction leaving. Ball intersect, `H` build and `Mwpm(H)` build
///                   are **discounted entirely**, exactly as `summarize_m2_artifact.py::critical_ns`
///                   discounts them: they are taken to run ahead of the syndrome they serve or to
///                   overlap the previous shot's matching. The certificate's terminal scan is *in*,
///                   because it decides whether the answer may be emitted at all and so cannot hide
///                   behind the next shot. That makes this an upper bound on an end-to-end win and an
///                   assumption about a machine that does not exist yet, which is why the end-to-end
///                   ratio (`stock_ns(G) / sum_total_ns`) is printed beside it and the stage split
///                   `intersect / h_build / mwpm_build / blossom_on_h / dual_scan / harvest` is
///                   reported in full.
///
///                   `q_this` rides along on every `(d, p, T)`, with `shots`, as the escalation rate
///                   the speedup campaign actually saw, and `escal_us = q_this * C_escalation` is the
///                   Phase-2 fallback amortised over all shots. §M3.2 caps what that can be worth at
///                   ~0.02% amortised, so it is charged beside the read rather than fitted: there is
///                   no separate deep single-point campaign, and no replay of the landed scheme's
///                   decision to difference against.
///
/// The error-rate grid runs **low**, `p ∈ {1e-4, 5e-4, 1e-3}`. That is where the ball graph is
/// meant to look best — fewer defects means a smaller `H` and a shorter solve against a stock decode
/// whose cost is set by the whole detector graph — so it is where the critical-path read has to be
/// taken. It also means every LER in the campaign is zero, which is exactly why §M7.6 level 3 is
/// read as a per-shot disagreement count rather than as a fit (§M6.4's vacuous-fit warning).
///
/// Usage:
///   two_phase_m7_artifact [--distances 5,7,9,11,13] [--error-rates 0.0001,0.0005,0.001]
///                         [--horizons 1.5,2.0] [--shots 20000] [--identity-shots 5000]
///                         [--mode scan|bitset] [--seed N] [--csv path]

#include <cstdio>

#include "benchmarks/two_phase/artifact_util.h"
#include "pyrematching/two_phase/driver/two_phase_decoding.h"

using namespace pm::two_phase;
using namespace pm::two_phase::artifact;

namespace {

struct Options {
    std::vector<size_t> distances = {5, 7, 9, 11, 13};
    std::vector<double> error_rates = {0.0001, 0.0005, 0.001};
    std::vector<double> horizons = {1.5, 2.0};
    size_t shots = 20000;
    size_t identity_shots = 5000;
    BallGraphBuildMode mode = BallGraphBuildMode::SCAN;
    uint64_t seed = 20260907;
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
        } else if (flag == "--identity-shots") {
            options.identity_shots = std::stoul(next());
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

TwoPhaseConfig config_for(double horizon_multiple, double unit, BallGraphBuildMode mode, bool stock_on_h) {
    TwoPhaseConfig config;
    config.T = horizon_multiple * unit;
    config.ball.T_max = horizon_multiple * unit;
    config.ball.R = 2.0 * config.ball.T_max;
    config.mode = mode;
    config.stock_on_h = stock_on_h;
    return config;
}

/// The `BallProfile` stage split, summed over a campaign. `TwoPhaseAggregateStats` does not carry
/// it — it is a level below — so it is read off `TwoPhaseDecoder::ball_profile`, which holds the
/// last decoded shot's split, once per shot.
struct StageSplit {
    long long intersect_ns{0};
    long long h_build_ns{0};
    long long mwpm_build_ns{0};
    long long blossom_on_h_ns{0};
    long long harvest_ns{0};
    long long dual_scan_ns{0};
    long long abandon_ns{0};
    long long total_ns{0};

    void add(const BallProfile& profile) {
        intersect_ns += profile.intersect_ns;
        h_build_ns += profile.h_build_ns;
        mwpm_build_ns += profile.mwpm_build_ns;
        blossom_on_h_ns += profile.blossom_on_h_ns;
        harvest_ns += profile.harvest_ns;
        dual_scan_ns += profile.dual_scan_ns;
        abandon_ns += profile.abandon_ns;
        total_ns += profile.total_ns;
    }
};

struct CampaignResult {
    TwoPhaseSummary summary;
    StageSplit split;
    long long sum_total_ns{0};

    /// §M2's critical path, applied to this front end: the stages that sit between a syndrome
    /// arriving and a correction leaving. `summarize_m2_artifact.py::critical_ns` defines it as
    /// blossom-on-`H` plus the readout, and takes ball intersect, `H` build, `Mwpm(H)` build and the
    /// unattributed remainder to be pipelined out — they either run ahead of the syndrome they serve
    /// or overlap the previous shot's matching. This is the denominator of the read.
    ///
    /// On §M7's front end the certificate joins them: the terminal scan sits between the solve and
    /// the emission and decides whether the answer may be emitted at all, so it cannot be hidden
    /// behind the next shot.
    ///
    /// `harvest_ns` is the readout on a certified shot and `abandon_shot` on an escalating one — the
    /// full charge either way, since both sit in the same place in the pipeline.
    inline long long crit_ns() const {
        return split.blossom_on_h_ns + split.dual_scan_ns + split.harvest_ns;
    }
};

/// One front end — stock blossom on `H` — and stock blossom on `G`, timed **on the same shot**.
///
/// `sum_stock_g_ns` is stock sparse blossom run to completion on `G`: no horizon, no ball graph, no
/// residual. It is the numerator of the read; `CampaignResult::crit_ns` is the denominator. Both
/// sides run the same solver, so nothing about the comparison depends on the truncated scheme —
/// which is not built here at all.
///
/// Per shot rather than through `decode_batch`, because the ball stage split has to be read off
/// `TwoPhaseDecoder::ball_profile` between shots.
///
/// The order within a shot is fixed — `H` then `G` — rather than rotated. A rotation is what an A/B
/// between two candidate front ends needs, and this is not one: the two sides here are measured on
/// different scales (a whole decode against a subset of one), a fixed position is the same for every
/// point in the grid, and an earlier three-way rotation was measured to *move* the operating point
/// rather than to centre it.
void run_campaign(
    const stim::DetectorErrorModel& dem,
    const TwoPhaseConfig& config,
    const std::vector<std::vector<uint64_t>>& shots,
    CampaignResult& out,
    long long& sum_stock_g_ns) {
    auto decoder = TwoPhaseDecoder::from_detector_error_model(dem, config, NUM_DISTINCT_WEIGHTS);
    // Stock exact decode's own instance, so that timing it does not disturb the front end's
    // `Mwpm(G)` — which the escalating path shares with the negative-weight preamble.
    pm::Mwpm reference = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
    std::vector<uint8_t> obs(std::max<size_t>(1, decoder.num_observables), 0);
    pm::total_weight_int weight = 0;

    // Warm both: the first shot is where the ball arenas, the `Mwpm(H)` node pool and stock's own
    // region arena reach their steady size, and it is not representative of any shot after it.
    if (!shots.empty()) {
        decoder.decode_to_obs(shots[0], obs.data(), weight, nullptr);
        pm::decode_detection_events(reference, shots[0], obs.data(), weight, false);
    }

    TwoPhaseProfile profile;
    HiResTimer timer;
    sum_stock_g_ns = 0;
    for (const auto& shot : shots) {
        decoder.decode_to_obs(shot, obs.data(), weight, &profile);
        decoder.stats.accumulate(profile);
        out.split.add(decoder.ball_profile);

        std::fill(obs.begin(), obs.end(), (uint8_t)0);
        weight = 0;
        timer.start();
        pm::decode_detection_events(reference, shot, obs.data(), weight, false);
        sum_stock_g_ns += timer.elapsed_ns();
    }
    out.summary = summarize(decoder.stats);
    out.sum_total_ns = decoder.stats.sum_total_ns;
}

/// §M7.6 level 1 + level 3, untimed: how many shots disagree with stock exact decode on `G` in
/// weight or in observable bytes. Expect 0. Kept out of the timed campaign entirely, because it
/// runs a second decoder on every shot.
struct IdentityResult {
    size_t shots{0};
    size_t weight_disagreements{0};
    size_t obs_disagreements{0};
    size_t certified{0};
    size_t escalated{0};
    size_t h_no_perfect_matching{0};
};

IdentityResult check_identity(
    const stim::DetectorErrorModel& dem,
    const TwoPhaseConfig& config,
    const std::vector<std::vector<uint64_t>>& shots) {
    auto decoder = TwoPhaseDecoder::from_detector_error_model(dem, config, NUM_DISTINCT_WEIGHTS);
    pm::Mwpm reference = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
    size_t num_observables = decoder.num_observables;

    IdentityResult result;
    std::vector<uint8_t> actual(std::max<size_t>(1, num_observables), 0);
    std::vector<uint8_t> expected(std::max<size_t>(1, num_observables), 0);
    for (const auto& shot : shots) {
        pm::total_weight_int actual_weight = 0;
        pm::total_weight_int expected_weight = 0;
        TwoPhaseProfile profile;
        decoder.decode_to_obs(shot, actual.data(), actual_weight, &profile);
        std::fill(expected.begin(), expected.end(), (uint8_t)0);
        pm::decode_detection_events(reference, shot, expected.data(), expected_weight, false);

        result.shots++;
        result.certified += (size_t)profile.certified;
        result.escalated += profile.escalated ? 1 : 0;
        result.h_no_perfect_matching += (size_t)profile.h_no_perfect_matching;
        if (actual_weight != expected_weight)
            result.weight_disagreements++;
        if (!std::equal(actual.begin(), actual.begin() + (long)num_observables, expected.begin()))
            result.obs_disagreements++;
    }
    return result;
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

    // ----------------------------------------------------------------------------- identity
    //
    // Read first and read as a count. A non-zero entry anywhere in this table voids everything
    // below it, so there is no point measuring speed until it is green.
    std::printf(
        "\n=== §M7.6 identity against stock exact decode on G (%zu shots per point, %s build) ===\n",
        options.identity_shots,
        mode_name);
    std::printf(
        "%4s %8s %5s %8s %10s %10s %10s %10s\n",
        "d",
        "p",
        "T",
        "shots",
        "wt_disagr",
        "obs_disagr",
        "certified",
        "no_pm");
    size_t total_disagreements = 0;
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, distance, noise, options.identity_shots, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);
            for (double horizon_multiple : options.horizons) {
                TwoPhaseConfig config = config_for(horizon_multiple, unit, options.mode, true);
                IdentityResult identity = check_identity(experiment.dem, config, experiment.shots);
                total_disagreements += identity.weight_disagreements + identity.obs_disagreements;
                std::printf(
                    "%4zu %8g %5g %8zu %10zu %10zu %10zu %10zu\n",
                    distance,
                    noise,
                    horizon_multiple,
                    identity.shots,
                    identity.weight_disagreements,
                    identity.obs_disagreements,
                    identity.certified,
                    identity.h_no_perfect_matching);
                emit("identity", distance, noise, horizon_multiple, "shots", (double)identity.shots);
                emit(
                    "identity",
                    distance,
                    noise,
                    horizon_multiple,
                    "weight_disagreements",
                    (double)identity.weight_disagreements);
                emit(
                    "identity",
                    distance,
                    noise,
                    horizon_multiple,
                    "obs_disagreements",
                    (double)identity.obs_disagreements);
                emit("identity", distance, noise, horizon_multiple, "certified", (double)identity.certified);
                emit("identity", distance, noise, horizon_multiple, "escalated", (double)identity.escalated);
                emit(
                    "identity",
                    distance,
                    noise,
                    horizon_multiple,
                    "h_no_perfect_matching",
                    (double)identity.h_no_perfect_matching);
            }
        }
    }
    std::printf("total disagreements across the grid: %zu (expected 0)\n", total_disagreements);

    // ----------------------------------------------------------------------------- speedup: the read
    //
    // `stock_us / crit_us`, with `crit_us = blossom_on_h + dual_scan + harvest`. Everything else the
    // shot does — ball intersect, `H` build, `Mwpm(H)` build, the unattributed remainder — is
    // discounted, exactly as `summarize_m2_artifact.py::critical_ns` discounts it. That is an
    // assumption about a machine that does not exist yet, not a measurement of this one, which is
    // why the end-to-end ratio is reported beside it and is the honest number for a CPU.
    std::printf(
        "\n=== the read: stock exact on G vs the M7 critical path on H, blossom + dual scan + harvest\n"
        "    (%zu shots per point, %s build) ===\n",
        options.shots,
        mode_name);
    std::printf(
        "%4s %8s %5s %10s %10s %10s %9s %9s %11s\n",
        "d",
        "p",
        "T",
        "stock_us",
        "crit_us",
        "e2e_us",
        "speedup",
        "e2e",
        "escal_us");
    std::vector<std::string> stage_rows;
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, distance, noise, options.shots, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double horizon_multiple : options.horizons) {
                CampaignResult stock;
                long long sum_stock_g_ns = 0;
                run_campaign(
                    experiment.dem,
                    config_for(horizon_multiple, unit, options.mode, true),
                    experiment.shots,
                    stock,
                    sum_stock_g_ns);

                double shots_run = (double)std::max<uint64_t>(1, stock.summary.shots);
                double stock_g_us = (double)sum_stock_g_ns / shots_run / 1000.0;
                double crit_us = (double)stock.crit_ns() / shots_run / 1000.0;
                double e2e_us = (double)stock.sum_total_ns / shots_run / 1000.0;
                // Escalation is a Phase-2 cost and sits outside the Phase-1 critical path the §M2
                // read is defined on, so it is amortised in separately rather than folded in:
                // `q * C_escalation` is what it adds to the mean shot.
                double escal_us = stock.summary.q * stock.summary.c_escalation / 1000.0;
                std::printf(
                    "%4zu %8g %5g %10.3f %10.3f %10.3f %9.2f %9.2f %11.4f\n",
                    distance,
                    noise,
                    horizon_multiple,
                    stock_g_us,
                    crit_us,
                    e2e_us,
                    crit_us > 0 ? stock_g_us / crit_us : 0.0,
                    e2e_us > 0 ? stock_g_us / e2e_us : 0.0,
                    escal_us);

                // The stage split, printed once the read is on the page. Percentages are of the
                // whole shot, so `blossom% + dscan% + hrvst%` is the share of it the read charges
                // and the rest is what the pipelining assumption discounts.
                double total = (double)std::max<long long>(1, stock.sum_total_ns);
                std::ostringstream row;
                row << distance << "," << noise << "," << horizon_multiple << ","
                    << 100.0 * (double)stock.split.intersect_ns / total << ","
                    << 100.0 * (double)stock.split.h_build_ns / total << ","
                    << 100.0 * (double)stock.split.mwpm_build_ns / total << ","
                    << 100.0 * (double)stock.split.blossom_on_h_ns / total << ","
                    << 100.0 * (double)stock.split.dual_scan_ns / total << ","
                    << 100.0 * (double)stock.split.harvest_ns / total;
                stage_rows.push_back(row.str());

                emit("speedup", distance, noise, horizon_multiple, "shots", shots_run);
                emit("speedup", distance, noise, horizon_multiple, "stock_g_ns", (double)sum_stock_g_ns);
                emit("speedup", distance, noise, horizon_multiple, "crit_ns", (double)stock.crit_ns());
                emit("speedup", distance, noise, horizon_multiple, "total_ns", (double)stock.sum_total_ns);
                emit("speedup", distance, noise, horizon_multiple, "mean_ns", stock.summary.measured_mean_ns);
                emit("speedup", distance, noise, horizon_multiple, "intersect_ns", (double)stock.split.intersect_ns);
                emit("speedup", distance, noise, horizon_multiple, "h_build_ns", (double)stock.split.h_build_ns);
                emit("speedup", distance, noise, horizon_multiple, "mwpm_build_ns", (double)stock.split.mwpm_build_ns);
                emit(
                    "speedup",
                    distance,
                    noise,
                    horizon_multiple,
                    "blossom_on_h_ns",
                    (double)stock.split.blossom_on_h_ns);
                emit("speedup", distance, noise, horizon_multiple, "dual_scan_ns", (double)stock.split.dual_scan_ns);
                emit("speedup", distance, noise, horizon_multiple, "harvest_ns", (double)stock.split.harvest_ns);
                emit("speedup", distance, noise, horizon_multiple, "q_this", stock.summary.q);
                emit(
                    "speedup",
                    distance,
                    noise,
                    horizon_multiple,
                    "amortised_escalation_ns",
                    stock.summary.q * stock.summary.c_escalation);
            }
        }
    }
    std::printf(
        "`speedup` is `stock_us / crit_us` and is an upper bound on the end-to-end win, not the win:\n"
        "it charges the shot only the stages assumed to sit between a syndrome arriving and a\n"
        "correction leaving, and compares against a stock decode that does the whole job in one place.\n"
        "`e2e` is `stock_us / e2e_us`, the same numerator over the whole two-phase shot, and is the\n"
        "honest number for a CPU. `escal_us` is `q_this * C_escalation`: the Phase-2 fallback\n"
        "amortised over all shots, which sits outside the Phase-1 critical path and is beside it, not\n"
        "in it. `q_this` itself is emitted per point with `shots`, so a zero there is this campaign's\n"
        "resolution floor rather than a measured absence of escalation.\n");

    std::printf("\n=== where the shot goes, as a share of the whole two-phase decode ===\n");
    std::printf(
        "%4s %8s %5s %9s %9s %9s %9s %9s %9s\n",
        "d",
        "p",
        "T",
        "isect%",
        "hbuild%",
        "mwpm%",
        "blossom%",
        "dscan%",
        "hrvst%");
    for (const std::string& row : stage_rows) {
        std::vector<double> values = parse_list<double>(row);
        if (values.size() != 9)
            continue;
        std::printf(
            "%4.0f %8g %5g %8.2f%% %8.2f%% %8.2f%% %8.2f%% %8.3f%% %8.2f%%\n",
            values[0],
            values[1],
            values[2],
            values[3],
            values[4],
            values[5],
            values[6],
            values[7],
            values[8]);
    }
    std::printf(
        "The last three columns are what the read charges; the first three are what it discounts as\n"
        "pipelined out, and are the whole distance between `speedup` and `e2e`.\n");

    std::printf("\ndone.\n");
    return total_disagreements == 0 ? 0 : 1;
}
