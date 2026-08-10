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

/// Per-shot latency of the two decoders, logged shot by shot rather than reduced to a mean.
///
/// The M7 artifact answers "how much faster, on average, at each grid point". This answers "what
/// does the latency *distribution* look like at one grid point", which is a different question and
/// needs a different output: one file per `(d, p, T, mode)`, one row per shot, and no aggregation at
/// all in the C++. `plot_latency_histograms.py` reads those files and draws the two histograms.
///
/// Two decoders, both **stock** sparse blossom, timed on the same shot in a fixed order:
///
///  - `stock_g_ns`  — stock exact decode on the original detector graph `G`. No horizon, no ball
///                    graph, no residual: the whole job in one place.
///  - the sparsified series — the §M7 front end on the sparsified graph `H`, **including
///                    escalation**. A shot the certificate rejects is charged its full Phase-2
///                    re-decode on `G`, so the escalation tail is in the histogram rather than
///                    amortised away beside it.
///
/// The sparsified series is not logged as one number: the three stages it is made of are logged
/// separately and summed by the reader.
///
///     sparse = blossom_ns + dscan_ns + hrvst_ns  (+ escal_stock_ns on an escalating shot)
///
/// which is the solve on `H`, §M7.7's terminal `max_u Y(u)` scan, and the harvest — the three stages
/// that are on the critical path once the front end is running — plus, on the shots that escalate,
/// the stock re-decode on `G` that Phase 2 pays. Ball intersect, `H` build and `Mwpm(H)` build are
/// **not** in it: §M2's critical-path read takes those three to run ahead of the syndrome they serve
/// or to overlap the previous shot's matching. They are not discarded either — their sum is logged
/// beside the rest as `excluded_ns`, so a reader who rejects the pipelining assumption can add them
/// back, and the plotting script's `--include-excluded` does exactly that.
///
/// Charging the stages rather than the shot's end-to-end `total_ns` means the unattributed cost
/// between stages is *not* in the series. That is the one thing this file's numbers are optimistic
/// about, and `total_ns` is logged as `total_ns` so the size of the gap stays checkable.
///
/// The order within a shot is fixed — `H` then `G` — rather than rotated, for the reason §M7.8
/// gives: the two sides are measured on different scales, a fixed position is the same at every
/// point in the grid, and a rotation was measured to move the operating point rather than centre it.
///
/// Shots the scheduler interfered with are flagged rather than dropped. `HiResTimer` is
/// thread-scoped, so such a shot is no longer *charged* for its off-CPU time, but it resumed on
/// caches it did not leave, so it is still an outlier — the plotting script excludes them by default
/// and says how many it excluded.
///
/// Usage:
///   two_phase_latency_profiler [--distances 5,7,9,11,13] [--error-rates 0.001]
///                              [--horizons 1.5,2.0] [--shots 20000] [--modes scan,bitset]
///                              [--warmup 256] [--seed N] [--out-dir DIR] [--tag NAME]
///
/// One log file per `(d, p, T, mode)`, named `latency_d{d}_p{p}_T{T}_{mode}.csv`, written to
/// `--out-dir`. See `write_header` for the schema.

#include <cstdio>
#include <filesystem>

#include "benchmarks/two_phase/artifact_util.h"
#include "pyrematching/two_phase/driver/two_phase_decoding.h"

using namespace pm::two_phase;
using namespace pm::two_phase::artifact;

namespace {

struct Options {
    std::vector<size_t> distances = {5, 7, 9, 11, 13};
    std::vector<double> error_rates = {0.001};
    std::vector<double> horizons = {1.5, 2.0};
    std::vector<BallGraphBuildMode> modes = {BallGraphBuildMode::SCAN};
    size_t shots = 20000;
    /// Shots decoded by both decoders before the timed loop starts, and not logged. The first shot
    /// is where the ball arenas, the `Mwpm(H)` node pool and stock's own region arena reach their
    /// steady size; a distribution read off a file whose first rows are allocation is not this
    /// decoder's distribution.
    size_t warmup = 256;
    uint64_t seed = 20260907;
    std::string out_dir = "benchmarks/two_phase/results/latency";
    std::string tag;
};

const char* mode_name(BallGraphBuildMode mode) {
    return mode == BallGraphBuildMode::BITSET ? "bitset" : "scan";
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
        } else if (flag == "--warmup") {
            options.warmup = std::stoul(next());
        } else if (flag == "--modes" || flag == "--mode") {
            options.modes.clear();
            for (const std::string& name : parse_list<std::string>(next())) {
                if (name == "bitset") {
                    options.modes.push_back(BallGraphBuildMode::BITSET);
                } else if (name == "scan") {
                    options.modes.push_back(BallGraphBuildMode::SCAN);
                } else {
                    throw std::invalid_argument("unrecognised isect mode " + name + " (want scan or bitset)");
                }
            }
            if (options.modes.empty())
                throw std::invalid_argument("--modes needs at least one of scan,bitset");
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--out-dir") {
            options.out_dir = next();
        } else if (flag == "--tag") {
            options.tag = next();
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
    // §M7's front end: stock blossom on `H` plus a max-dual certificate. The other side of the
    // comparison is the same solver on `G`, so nothing here depends on the truncated scheme.
    config.stock_on_h = true;
    // §M6.4's timer discipline. The flag costs one `getrusage` per shot and is what lets a
    // descheduled shot be flagged instead of silently widening the tail.
    config.detect_preemption = true;
    return config;
}

/// One shot's row. The sparsified side is stored stage by stage rather than pre-summed, so the
/// reader can see what it is made of and the definition lives in one visible place.
struct ShotRow {
    long long stock_g_ns{0};
    /// The solve on `H`. Excludes the dual scan, so the two are additive.
    long long blossom_ns{0};
    /// §M7.7's terminal `max_u Y(u)` scan: the only cost the certificate adds.
    long long dscan_ns{0};
    long long hrvst_ns{0};
    /// The Phase-2 stock re-decode on `G`, on escalating shots only; 0 otherwise. Timed inside the
    /// escalation call, so it is that shot's own re-decode rather than the mean of the column beside
    /// it.
    long long escal_stock_ns{0};
    /// intersect + `H` build + `Mwpm(H)` build: the stages discounted as pipelined out.
    long long excluded_ns{0};
    /// The whole two-phase shot end to end. Not charged anywhere; kept so that the gap between it
    /// and the stages above stays visible.
    long long total_ns{0};
    int defects{0};
    int escalated{0};
    int certified{0};
    int contaminated{0};

    /// The sparsified series. Kept identical to `plot_latency_histograms.py`'s `series_of`.
    long long sparse_ns() const {
        return blossom_ns + dscan_ns + hrvst_ns + escal_stock_ns;
    }
};

std::string point_name(size_t distance, double noise, double horizon, BallGraphBuildMode mode, const std::string& tag) {
    char buffer[256];
    std::snprintf(
        buffer,
        sizeof(buffer),
        "latency_d%zu_p%g_T%g_%s%s%s.csv",
        distance,
        noise,
        horizon,
        mode_name(mode),
        tag.empty() ? "" : "_",
        tag.c_str());
    return buffer;
}

/// The schema, and everything needed to read a row without the command line that produced it.
///
/// `#` lines are `key=value` metadata; the first non-`#` line is the column header. The timer
/// backend and whether it is thread-scoped are recorded because a campaign that fell back to a
/// wall-clock backend must not be read as if it had not (§M6.1).
void write_header(
    std::ofstream& out,
    const Options& options,
    size_t distance,
    double noise,
    double horizon,
    BallGraphBuildMode mode,
    double unit) {
    out << "# schema=two_phase_latency_v2\n";
    out << "# d=" << distance << "\n";
    out << "# p=" << noise << "\n";
    out << "# T=" << horizon << "\n";
    out << "# T_weight_units=" << horizon * unit << "\n";
    out << "# mode=" << mode_name(mode) << "\n";
    out << "# shots=" << options.shots << "\n";
    out << "# warmup=" << options.warmup << "\n";
    out << "# seed=" << options.seed << "\n";
    out << "# rounds=" << distance << "\n";
    out << "# circuit=rotated_memory_x\n";
    out << "# timer_backend=" << pm::perf::current_backend_name() << "\n";
    out << "# timer_thread_scoped=" << (pm::perf::backend_is_thread_scoped(pm::perf::current_backend()) ? 1 : 0)
        << "\n";
    // What each latency column is, stated in the file rather than left to the plotting script.
    out << "# stock_g_ns=stock exact sparse blossom on the original detector graph G\n";
    out << "# sparse_ns=blossom_ns+dscan_ns+hrvst_ns+escal_stock_ns, summed by the reader\n";
    out << "# blossom_ns=the solve on the sparsified graph H, dual scan excluded\n";
    out << "# dscan_ns=the terminal max_u Y(u) scan, i.e. the certificate's own cost\n";
    out << "# hrvst_ns=harvest/extraction on H\n";
    out << "# escal_stock_ns=Phase-2 stock re-decode on G, on escalating shots only, else 0\n";
    out << "# excluded_ns=intersect+h_build+mwpm_build, the stages discounted as pipelined out\n";
    out << "# total_ns=the whole two-phase shot end to end; not charged, kept so the gap to the"
           " stage sum stays visible\n";
    out << "shot,defects,escalated,certified,contaminated,stock_g_ns,blossom_ns,dscan_ns,hrvst_ns,"
           "escal_stock_ns,excluded_ns,total_ns\n";
}

/// Both decoders on every shot, `H` first. Returns one row per shot, in shot order.
std::vector<ShotRow> run_point(
    const stim::DetectorErrorModel& dem,
    const TwoPhaseConfig& config,
    const std::vector<std::vector<uint64_t>>& shots,
    size_t warmup) {
    auto decoder = TwoPhaseDecoder::from_detector_error_model(dem, config, NUM_DISTINCT_WEIGHTS);
    // Stock exact decode gets its own instance, so timing it does not disturb the front end's
    // `Mwpm(G)` — which the escalating path shares with the negative-weight preamble.
    pm::Mwpm reference = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
    std::vector<uint8_t> obs(std::max<size_t>(1, decoder.num_observables), 0);
    pm::total_weight_int weight = 0;

    for (size_t i = 0; i < std::min(warmup, shots.size()); i++) {
        decoder.decode_to_obs(shots[i], obs.data(), weight, nullptr);
        pm::decode_detection_events(reference, shots[i], obs.data(), weight, false);
    }

    std::vector<ShotRow> rows;
    rows.reserve(shots.size());
    TwoPhaseProfile profile;
    HiResTimer timer;
    for (const auto& shot : shots) {
        ShotRow row;

        decoder.decode_to_obs(shot, obs.data(), weight, &profile);
        const BallProfile& ball = decoder.ball_profile;
        row.excluded_ns = ball.intersect_ns + ball.h_build_ns + ball.mwpm_build_ns;
        row.blossom_ns = ball.blossom_on_h_ns;
        row.dscan_ns = ball.dual_scan_ns;
        row.hrvst_ns = ball.harvest_ns;
        // Zero unless the shot escalated: `stock_ns` is written by `escalate_to_stock` and by
        // nothing else.
        row.escal_stock_ns = profile.stock_ns;
        row.total_ns = profile.total_ns;
        row.defects = profile.num_defects;
        row.escalated = profile.escalated ? 1 : 0;
        row.certified = profile.certified;
        row.contaminated = profile.contaminated ? 1 : 0;

        // Stock on `G`, on the same shot, with its own preemption probe: either side being
        // descheduled makes the pair an outlier, so the flag covers both.
        PreemptionProbe before;
        before.sample();
        std::fill(obs.begin(), obs.end(), (uint8_t)0);
        weight = 0;
        timer.start();
        pm::decode_detection_events(reference, shot, obs.data(), weight, false);
        row.stock_g_ns = timer.elapsed_ns();
        PreemptionProbe after;
        after.sample();
        if (after.switched_since(before))
            row.contaminated = 1;

        rows.push_back(row);
    }
    return rows;
}

/// The mean over the uncontaminated rows, printed to stdout so a run says what it wrote without
/// waiting for the plots. The files are the artifact; this is a progress line.
///
/// The mean is taken over *every* uncontaminated shot, escalated ones included — an escalated shot
/// is charged its Phase-2 re-decode on `G` inside the sparsified series, so dropping it would price
/// the front end at a rate no deployment gets. That is what makes the printed ratio a speedup rather
/// than a best-case.
struct QuickStats {
    double mean{0};
    size_t kept{0};
};

template <typename Charge>
QuickStats quick_stats(const std::vector<ShotRow>& rows, Charge charge) {
    QuickStats stats;
    double total = 0;
    for (const ShotRow& row : rows) {
        if (row.contaminated)
            continue;
        total += (double)charge(row);
        stats.kept++;
    }
    if (stats.kept)
        stats.mean = total / (double)stats.kept;
    return stats;
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

    std::printf(
        "timer backend: %s (%s)\n",
        pm::perf::current_backend_name(),
        pm::perf::backend_is_thread_scoped(pm::perf::current_backend()) ? "thread-scoped"
                                                                        : "WALL CLOCK — timings include off-CPU time");
    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }
    std::printf("writing one log file per (d, p, T, mode) to %s\n\n", options.out_dir.c_str());
    std::printf(
        "%4s %8s %5s %8s %8s %11s %11s %8s %8s %9s %8s\n",
        "d",
        "p",
        "T",
        "mode",
        "shots",
        "stock_mean",
        "sparse_mean",
        "speedup",
        "escal",
        "escal_frac",
        "contam");

    size_t files_written = 0;
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            Experiment experiment = generate(distance, distance, noise, options.shots, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double horizon : options.horizons) {
                for (BallGraphBuildMode mode : options.modes) {
                    std::vector<ShotRow> rows =
                        run_point(experiment.dem, config_for(horizon, unit, mode), experiment.shots, options.warmup);

                    std::string path = options.out_dir + "/" + point_name(distance, noise, horizon, mode, options.tag);
                    std::ofstream out(path);
                    if (!out.is_open()) {
                        std::cerr << "error: could not open " << path << " for writing"
                                  << " (does " << options.out_dir << " exist?)\n";
                        return 1;
                    }
                    write_header(out, options, distance, noise, horizon, mode, unit);
                    for (size_t i = 0; i < rows.size(); i++) {
                        const ShotRow& row = rows[i];
                        out << i << "," << row.defects << "," << row.escalated << "," << row.certified << ","
                            << row.contaminated << "," << row.stock_g_ns << "," << row.blossom_ns << ","
                            << row.dscan_ns << "," << row.hrvst_ns << "," << row.escal_stock_ns << ","
                            << row.excluded_ns << "," << row.total_ns << "\n";
                    }
                    out.close();
                    files_written++;

                    QuickStats stock = quick_stats(rows, [](const ShotRow& row) { return row.stock_g_ns; });
                    QuickStats sparse = quick_stats(rows, [](const ShotRow& row) { return row.sparse_ns(); });
                    size_t escalated = 0;
                    size_t contaminated = 0;
                    for (const ShotRow& row : rows) {
                        escalated += (size_t)row.escalated;
                        contaminated += (size_t)row.contaminated;
                    }
                    // stock's mean over the two-phase mean, both including every escalated shot:
                    // what the front end buys once the escalation tail is paid for.
                    double speedup = sparse.mean > 0 ? stock.mean / sparse.mean : 0.0;
                    double escalated_fraction = rows.empty() ? 0.0 : (double)escalated / (double)rows.size();
                    std::printf(
                        "%4zu %8g %5g %8s %8zu %9.3fus %9.3fus %8.3fx %8zu %8.3f%% %8zu\n",
                        distance,
                        noise,
                        horizon,
                        mode_name(mode),
                        sparse.kept,
                        stock.mean / 1000.0,
                        sparse.mean / 1000.0,
                        speedup,
                        escalated,
                        100.0 * escalated_fraction,
                        contaminated);
                    std::fflush(stdout);
                }
            }
        }
    }

    std::printf(
        "\n%zu file(s) written. Means above are over uncontaminated shots only, escalated shots"
        " included;\n`contam` is how many of the %zu shots at each point were dropped from them, and"
        " `escal`/`escal_frac`\nare over all %zu. `speedup` is stock_mean/sparse_mean, so the"
        " escalation tail is priced into it.\n`sparse_mean` is blossom+dscan+hrvst, plus the Phase-2"
        " re-decode on the shots that escalate; ball\nintersect, H build and Mwpm(H) build are not in"
        " it and are logged as `excluded_ns` in every file.\n",
        files_written,
        options.shots,
        options.shots);
    return 0;
}
