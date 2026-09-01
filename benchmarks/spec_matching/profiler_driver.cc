// Copyright 2026 SpecMatching contributors
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
/// A mean answers "how much faster, on average, at each grid point". This answers "what does the
/// latency *distribution* look like at one grid point", which is a different question and needs a
/// different output: one file per `(d, p, T, mode)`, one row per shot, and no aggregation at all in
/// the C++. `plot_latency_histograms.py` reads those files and draws the two histograms.
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
/// ## The two things the report prices, and which one the speedup is
///
/// `sparse` above is one **serial** machine: it runs the front end on `H`, and on a shot the
/// certificate rejects it *then* pays a stock re-decode on `G`. That is the honest price of a
/// decoder that only starts `G` once it knows it needs it, and it is reported as `sparse_mean`.
///
/// It is not the machine anyone would deploy. A deployment has both graphs in front of it and no
/// reason to hold `G` back, so it starts the solve on `H` and the solve on `G` **concurrently** and
/// stops the moment either produces a matching it can use. That system is what `system_mean` prices:
///
///     system = stock_g_ns                                  on an escalating shot
///            = min(blossom_ns + dscan_ns + hrvst_ns,
///                  stock_g_ns)                             otherwise
///
/// On a shot the certificate keeps, the answer is whichever race finished first — in practice `H`'s,
/// which is the whole point, but the `min` is written rather than assumed. On a shot it rejects,
/// `H`'s work bought nothing and the system waits out the `G` decode that has been running all
/// along: the escalation costs the *full* `G` latency and not a nanosecond more, where the serial
/// machine pays `H` and then `G` end to end. Escalation is therefore priced as lost opportunity
/// rather than as an added tail, which is what concurrency actually buys.
///
/// Every speedup in this report is `stock_mean / system_mean` — the decoder system without graph
/// sparsification against the same system predicating on `H` — over the same uncontaminated shots.
/// `sparse_mean` is printed beside it as the serial reading, and is never a denominator.
///
/// Two costs the concurrent reading assumes away, stated rather than buried: the two solves are
/// assumed to run on cores that are not competing for each other's memory bandwidth, and `G`'s
/// decode is assumed to be running from the start of the shot at no scheduling cost. Both are the
/// same critical-path style of assumption §M2 makes for ball intersect and the `H` build.
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
/// caches it did not leave, so it is still an outlier — the printed means exclude them, and the
/// plotting script excludes them by default and says how many it excluded.
///
/// Both the timer and that flag are platform-dependent, and the run says which instruments it got
/// before it decodes anything. On Linux the timer is a per-thread PMU cycle counter and the flag is
/// the thread's context-switch counters; on macOS both come from
/// `clock_gettime(CLOCK_THREAD_CPUTIME_ID)` — the timer directly, the flag as the gap between wall
/// time and thread CPU time over the shot, since Darwin exposes no per-thread switch counter (see
/// `PreemptionProbe`). If neither a thread-scoped clock is available, the run prints a wall-clock
/// warning banner and every number it produces is an upper bound with the scheduler inside it.
///
/// Beside the latency series, every run also dumps the **component structure** of the sparsified
/// graph `H`: its connected components, their sizes, weighted diameters and boundary structure, the
/// `H` edge-weight and boundary-cost distributions, and how much of the defect set sits in a
/// component of size `<= k`, resolved without the solver at all. That is a structural
/// measurement, not a timing one — it is computed after each shot's timed window has closed and is
/// charged to nothing — and it goes to its own file per grid point, so the latency schema keeps its
/// meaning and the two are read together.
///
/// `blossom_ns` and `hrvst_ns` are the latency figures of merit to read beside it: the solve on `H`
/// and the harvest are what a smaller node set moves, so they are summed into `b+h` in the printed
/// table as well as logged per shot.
///
/// `--k` takes a list — `--k 0,1,2,3,4`, or `--k 3,1,2,4` — and runs one **independent experiment**
/// per `k` in it, in the order given, on the same shots. `k` is the upper bound on the size of a
/// connected component of `H` the small-component resolver takes off the solver: the solver sees
/// only components of size `> k`, and the sizes at or below it are resolved exactly, before it, in
/// series. That resolve is a real cost on the critical path and it is **not** in any latency column
/// here: this driver measures solver + harvest on the size-`> k` graph, and the resolve is outside
/// that scope by intent, to be measured separately. `k` is capped at 4, where the resolver stops
/// being defined.
///
/// Usage:
///   profiler_driver [--distances 5,7,9,11,13] [--error-rates 0.001]
///                              [--horizons 1.5,2.0] [--k 0,1,2,3,4] [--shots 20000]
///                              [--modes scan,bitset] [--warmup 256] [--seed N] [--out-dir DIR]
///                              [--tag NAME] [--progress auto|always|never] [--no-component-stats]
///                              [--no-prune]
///
/// One log file per `(d, p, T, mode, k)`, named `latency_d{d}_p{p}_T{T}_{mode}_k{k}.csv`, written to
/// `--out-dir`, and one component file beside it named `components_...csv`. See `write_header` and
/// `write_component_file` for the two schemas.
///
/// A progress bar is drawn per grid point when stderr is a terminal (`ProgressBar`, `--progress`).
/// It says the run is alive; it does not keep the machine awake. On a laptop, hold the wake lock
/// separately — `caffeinate -dimsu profiler_driver ...` on macOS — because a sweep whose
/// second half is measured in a different power state than its first is two distributions.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#include "benchmarks/spec_matching/profiler_util.h"
#include "specmatching/spec_matching/driver/spec_matching_decoding.h"

using namespace pm::spec_matching;
using namespace pm::spec_matching::profiler;

namespace {

enum class ProgressMode { AUTO, ALWAYS, NEVER };

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
    std::string out_dir = "benchmarks/spec_matching/results/latency";
    std::string tag;
    /// `AUTO` draws the bar when stderr is a terminal and stays quiet when it is not, so a
    /// redirected run's log is not a megabyte of carriage returns. `always` forces it on, and into
    /// a file or a pipe it degrades to one plain line per 10% — which is what a `nohup`'d overnight
    /// run wants.
    ProgressMode progress = ProgressMode::AUTO;
    /// The component structure of `H`, per shot. On by default: it is a first-class output of this
    /// driver, not a diagnostic. `--no-component-stats` turns it off so the same binary can produce
    /// a latency-only run to compare against.
    bool component_stats = true;
    /// §A's `k`, one **independent experiment** per entry, run in the order given on the same shots.
    /// The solver sees only components of size `> k`; sizes `1..k` are resolved exactly off it by
    /// the small-component resolver, in series, before it — a real cost that is deliberately in no
    /// latency column here (see the file comment). The figures of merit are `blossom_ns` and
    /// `hrvst_ns` — the solver and the harvest on the size-`> k` graph — and their aggregates.
    ///
    /// `k = 0` resolves nothing and is the un-pruned path, which is what the rest are read against;
    /// `k = 2` is the previous branch's production behaviour. Capped at 4, where the resolver stops
    /// being defined. `--no-prune` is shorthand for `--k 0`.
    std::vector<int> ks = {2};
};

/// The largest `k` a resolver is defined for; the decoder rejects anything above it at
/// construction, and this driver rejects it at the flag so the error names the flag.
const int MAX_K = 4;

const char* mode_name(BallGraphBuildMode mode) {
    return mode == BallGraphBuildMode::BITSET ? "bitset" : "scan";
}

bool stderr_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    if (isatty(fileno(stderr)) != 1)
        return false;
    const char* term = std::getenv("TERM");
    // A dumb terminal takes neither the carriage return nor the erase-to-end-of-line below.
    return term != nullptr && std::string(term) != "dumb";
#else
    return false;
#endif
}

std::string format_duration(double seconds) {
    if (!(seconds >= 0) || seconds > 1e7)
        return "--";
    char buffer[32];
    int whole = (int)(seconds + 0.5);
    if (whole >= 3600)
        std::snprintf(buffer, sizeof(buffer), "%dh%02dm", whole / 3600, (whole % 3600) / 60);
    else if (whole >= 60)
        std::snprintf(buffer, sizeof(buffer), "%dm%02ds", whole / 60, whole % 60);
    else
        std::snprintf(buffer, sizeof(buffer), "%ds", whole);
    return buffer;
}

/// A per-job progress bar on stderr, sized so that drawing it cannot show up in the numbers.
///
/// What it is for: a sweep left running unattended on a laptop is a sweep whose second half may be
/// measured on a machine in a different power state than its first half, and a distribution taken
/// across that boundary is two distributions. A drawn bar is not itself a wake lock — on macOS,
/// `caffeinate -dimsu profiler_driver ...` is what actually holds the machine up — but it
/// is what tells a run that is still moving apart from one that has stalled, and it keeps the
/// terminal from being a dark, idle window for an hour.
///
/// It cannot perturb what is being measured, by construction:
///
///  - `tick` is called between shots only, after the stock decode's second `PreemptionProbe` sample
///    and before the next shot's decode takes its first. Every timer is stopped and every probe
///    interval is closed when it runs, so a redraw that blocks on the terminal lands in no
///    measurement window and cannot flag a shot as contaminated either.
///  - On all but one shot in `total/100` (one in `total/10` when not a terminal), `tick` is a
///    compare against a precomputed threshold and a return. The clock reading and the write happen
///    at most ~101 times per job, against 20000-odd shots.
///  - It writes to stderr, so the summary table on stdout survives a redirect unmixed, and the bar
///    is still on the screen when stdout is redirected away.
struct ProgressBar {
    static const size_t WIDTH = 32;

    bool enabled{false};
    bool interactive{false};
    std::string label;
    size_t total{0};
    /// Shots between redraws. Precomputed so the per-shot cost is one comparison.
    size_t step{1};
    size_t next_redraw{0};
    std::chrono::steady_clock::time_point started;

    /// Starts a job's bar. `label` — the job's grid point — is set by the caller beforehand.
    void begin(size_t job_total) {
        total = job_total;
        step = std::max<size_t>(1, job_total / (interactive ? 100 : 10));
        next_redraw = step;
        started = std::chrono::steady_clock::now();
        if (enabled && interactive)
            draw(0);
    }

    /// Called once per shot from inside the loop; `done` counts shots retired, warmup included.
    void tick(size_t done) {
        if (!enabled || done < next_redraw)
            return;
        next_redraw = done + step;
        draw(done);
    }

    /// Interactively, the line is erased and left to whatever is printed next: the job's row on
    /// stdout is the record, and a half-drawn bar above it is not worth keeping. In a log, the job
    /// gets one closing line at 100% so its wall time is on the record too.
    void finish() {
        if (!enabled)
            return;
        if (interactive)
            std::fprintf(stderr, "\r\x1b[K");
        else
            draw(total);
        std::fflush(stderr);
    }

    /// A line that is not part of a bar — what the run is doing while no job is in flight.
    void note(const std::string& text) {
        if (!enabled)
            return;
        std::fprintf(stderr, "%s%s\n", interactive ? "\r\x1b[K" : "", text.c_str());
        std::fflush(stderr);
    }

    void draw(size_t done) {
        double fraction = total > 0 ? (double)done / (double)total : 1.0;
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
        double left = fraction > 0.0 ? elapsed * (1.0 - fraction) / fraction : -1.0;
        if (!interactive) {
            std::fprintf(
                stderr,
                "%s  %3.0f%%  %zu/%zu  %s elapsed\n",
                label.c_str(),
                100.0 * fraction,
                done,
                total,
                format_duration(elapsed).c_str());
            std::fflush(stderr);
            return;
        }
        char cells[WIDTH + 1];
        size_t filled = std::min(WIDTH, (size_t)(fraction * (double)WIDTH));
        for (size_t i = 0; i < WIDTH; i++)
            cells[i] = i < filled ? '#' : '.';
        cells[WIDTH] = '\0';
        // `\x1b[K` erases the tail of the previous, possibly longer, line; `\r` alone would leave it.
        std::fprintf(
            stderr,
            "\r%s [%s] %3.0f%%  %zu/%zu  %s elapsed, %s left\x1b[K",
            label.c_str(),
            cells,
            100.0 * fraction,
            done,
            total,
            format_duration(elapsed).c_str(),
            format_duration(left).c_str());
        std::fflush(stderr);
    }
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
        } else if (flag == "--progress") {
            std::string value = next();
            if (value == "auto") {
                options.progress = ProgressMode::AUTO;
            } else if (value == "always" || value == "on") {
                options.progress = ProgressMode::ALWAYS;
            } else if (value == "never" || value == "off") {
                options.progress = ProgressMode::NEVER;
            } else {
                throw std::invalid_argument("unrecognised --progress " + value + " (want auto, always or never)");
            }
        } else if (flag == "--no-progress") {
            options.progress = ProgressMode::NEVER;
        } else if (flag == "--no-component-stats") {
            options.component_stats = false;
        } else if (flag == "--k") {
            // The list's order is the run order: `--k 3,1,2,4` runs them in that order, and each
            // one is an independent experiment over the same shots. Duplicates are rejected rather
            // than deduplicated, because two runs at one `k` would write the same file twice and
            // the second would silently be the one kept.
            options.ks = parse_list<int>(next());
            if (options.ks.empty())
                throw std::invalid_argument("--k needs at least one value");
            for (size_t a = 0; a < options.ks.size(); a++) {
                int k = options.ks[a];
                if (k < 0 || k > MAX_K)
                    throw std::invalid_argument(
                        "--k value " + std::to_string(k) + " is outside [0, " + std::to_string(MAX_K) +
                        "]; no small-component resolver is defined above " + std::to_string(MAX_K) + ".");
                for (size_t b = 0; b < a; b++) {
                    if (options.ks[b] == k)
                        throw std::invalid_argument("--k lists " + std::to_string(k) + " twice");
                }
            }
        } else if (flag == "--no-prune") {
            options.ks = {0};
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    return options;
}

SpecMatchingConfig config_for(
    double horizon_multiple, double unit, BallGraphBuildMode mode, bool component_stats, int k) {
    SpecMatchingConfig config;
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
    // §C. The decomposition runs after each shot's `total_ns` has been read and is timed by
    // nothing, so it cannot enter any column of the latency file — which is why the standard run
    // can leave it on rather than needing a second pass.
    config.collect_component_stats = component_stats;
    // §A. The resolve is a serial pre-pass on the critical path and is charged to no column of the
    // latency file — by measurement scope, not because it is free: what this driver reports is
    // `blossom_on_h_ns` and `harvest_ns` on the size-`> k` graph.
    config.prune_component_max_size = k;
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
    /// The whole spec-matching shot end to end. Not charged anywhere; kept so that the gap between it
    /// and the stages above stays visible.
    long long total_ns{0};
    int defects{0};
    int escalated{0};
    int certified{0};
    int contaminated{0};

    /// §C.1, on the shots where it was collected. Structural, not timed: these describe the graph
    /// the stages above ran on, and none of them is in any latency column.
    ///
    /// `h_defects` is `H`'s node count — the post-preamble defects, which is what the component
    /// counts partition and therefore the denominator any per-shot fraction of them is taken over.
    /// It is *not* `defects`, which is the raw detection-event count.
    int h_defects{0};
    int components{0};
    int singleton_components{0};
    int pair_components{0};
    int nontrivial_components{0};
    int largest_component{0};
    int defects_in_trivial{0};
    int defects_to_solver{0};

    /// §A's per-shot run label: `H` defects the small-component resolver settled off the solver at
    /// this run's `k`, sizes `1..k` included. Structural, not timed — the resolve is a serial
    /// pre-pass excluded from every column above by measurement scope — and 0 at `k = 0`.
    int defects_resolved_small{0};

    /// The sparsified series, as one **serial** machine: the front end on `H`, then — on a shot the
    /// certificate rejects — the Phase-2 re-decode on `G` after it. Kept identical to
    /// `plot_latency_histograms.py`'s `series_of`.
    long long sparse_ns() const {
        return blossom_ns + dscan_ns + hrvst_ns + escal_stock_ns;
    }

    /// The front end alone: what `H` costs before anything is known about whether it will be kept.
    long long h_path_ns() const {
        return blossom_ns + dscan_ns + hrvst_ns;
    }

    /// The **concurrent** decoder system: `H` and `G` started together, the shot ending the moment a
    /// usable matching exists. See the file comment for the definition and its assumptions.
    ///
    /// On an escalating shot that is the `G` decode, which has been running all along — `escal_stock_ns`
    /// is deliberately *not* added, because in this system the re-decode is not a second decode. On
    /// any other shot it is whichever of the two finished first; the `min` is written out rather than
    /// assumed, so a grid point where `H` is the slower side reports that instead of hiding it.
    long long system_ns() const {
        if (escalated)
            return stock_g_ns;
        return std::min(h_path_ns(), stock_g_ns);
    }
    /// §D's figure of merit for the pruning experiment: the two stages a smaller node set moves.
    long long solve_and_harvest_ns() const {
        return blossom_ns + hrvst_ns;
    }
};

/// `k` is in the name, not only in the header: the runs of a `--k` list differ in nothing else, so
/// without it the second experiment would overwrite the first.
std::string point_name(
    const char* prefix,
    size_t distance,
    double noise,
    double horizon,
    BallGraphBuildMode mode,
    int k,
    const std::string& tag) {
    char buffer[256];
    std::snprintf(
        buffer,
        sizeof(buffer),
        "%s_d%zu_p%g_T%g_%s_k%d%s%s.csv",
        prefix,
        distance,
        noise,
        horizon,
        mode_name(mode),
        k,
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
    int k,
    double unit,
    const ProgressBar& progress) {
    out << "# schema=spec_matching_latency_v4\n";
    out << "# d=" << distance << "\n";
    out << "# p=" << noise << "\n";
    out << "# T=" << horizon << "\n";
    out << "# k=" << k << "\n";
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
    // How `contaminated` was decided on the machine that wrote this file. A rate of zero means two
    // different things under `rusage_thread_switches` and under `none`, and only the header can tell
    // them apart (§M6.4).
    out << "# preemption_probe=" << PreemptionProbe::mechanism_name() << "\n";
    out << "# preemption_slack_ns=" << PreemptionProbe::slack_ns() << "\n";
    // Recorded for the same reason the backend is: it is a write this process made while the loop
    // was running. It happens between shots, outside every timed window, at most ~101 times a job.
    out << "# progress_bar=" << (progress.enabled ? (progress.interactive ? "bar" : "lines") : "off") << "\n";
    // What each latency column is, stated in the file rather than left to the plotting script.
    out << "# stock_g_ns=stock exact sparse blossom on the original detector graph G\n";
    out << "# sparse_ns=blossom_ns+dscan_ns+hrvst_ns+escal_stock_ns, summed by the reader; the"
           " serial machine, which starts G only once H is rejected\n";
    // The system the report's speedup is a ratio of. Derived by the reader from columns that are all
    // already here, so it is a definition rather than a column — stated in the file so that every
    // reader derives the same one.
    out << "# system_ns=the concurrent decoder system: H and G solved at once, the shot ending when a"
           " usable matching exists.\n";
    out << "#   escalated shot: stock_g_ns (the G decode was already running; H bought nothing)\n";
    out << "#   otherwise:      min(blossom_ns+dscan_ns+hrvst_ns, stock_g_ns)\n";
    out << "# speedup=mean(stock_g_ns)/mean(system_ns) over uncontaminated shots\n";
    out << "# blossom_ns=the solve on the sparsified graph H, dual scan excluded\n";
    out << "# dscan_ns=the terminal max_u Y(u) scan, i.e. the certificate's own cost\n";
    out << "# hrvst_ns=harvest/extraction on H\n";
    out << "# escal_stock_ns=Phase-2 stock re-decode on G, on escalating shots only, else 0\n";
    out << "# excluded_ns=intersect+h_build+mwpm_build, the stages discounted as pipelined out\n";
    out << "# total_ns=the whole spec-matching shot end to end; not charged, kept so the gap to the"
           " stage sum stays visible\n";
    // §D's figure of merit, named in the file so a reader does not have to be told which columns to
    // add: the prune moves the solver and the harvest, and nothing else in this schema.
    out << "# figure_of_merit=blossom_ns+hrvst_ns\n";
    // Which experiment of the `--k` list this file is. The runs differ in `k` and in nothing else,
    // and each is independent; what moves between them is the figure of merit above.
    out << "# prune_component_max_size=" << k << "\n";
    out << "# k_semantics=the solver sees only components of size > k; sizes 1..k are resolved"
           " exactly off it\n";
    // Stated in the file, because a reader who adds the resolve back has to know it was never in
    // here — and that leaving it out is a scope decision about what this branch measures, not a
    // claim that the resolve is free or concurrent.
    out << "# resolve_excluded=the size <= k resolve is a serial pre-pass on the critical path,"
           " excluded from every latency column here by measurement scope\n";
    out << "# defects_resolved_small=H defects the resolver settled off the solver this shot"
           " (sizes 1..k); structural, untimed\n";
    // The component columns of schema v3. Structural and untimed — they describe `H`, they are not
    // part of any series, and they are all 0 under `--no-component-stats`.
    out << "# component_stats=" << (options.component_stats ? 1 : 0) << "\n";
    out << "# h_defects=nodes of H, i.e. the post-preamble defects the components partition\n";
    out << "# components=connected components of H over its defect-defect edges\n";
    out << "# defects_in_trivial=defects in components of size <= 2\n";
    out << "# defects_to_solver=defects a size <= 2 resolver would have had to leave to the solver\n";
    out << "shot,defects,escalated,certified,contaminated,stock_g_ns,blossom_ns,dscan_ns,hrvst_ns,"
           "escal_stock_ns,excluded_ns,total_ns,"
           "h_defects,components,singleton_components,pair_components,nontrivial_components,"
           "largest_component,defects_in_trivial,defects_to_solver,defects_resolved_small\n";
}

/// §C.2/§C.3's distributions for one grid point, in long form so that one schema covers scalars and
/// histograms alike: `kind` names the table, `key` the statistic or the bin, `value` the count.
///
/// Bins of the three weight histograms are in sixteenths of `T` — bin `k` is
/// `[k * T / bins_per_T, (k + 1) * T / bins_per_T)`, last bin overflowing — which is recorded in the
/// metadata rather than left to the reader to know. Sizes and degrees are binned by count, with the
/// last bin overflowing the same way.
void write_component_file(
    std::ofstream& out,
    const Options& options,
    size_t distance,
    double noise,
    double horizon,
    BallGraphBuildMode mode,
    int k,
    double unit,
    horizon_int horizon_time_units,
    const SpecMatchingSummary& summary) {
    out << "# schema=spec_matching_components_v1\n";
    out << "# d=" << distance << "\n";
    out << "# p=" << noise << "\n";
    out << "# T=" << horizon << "\n";
    out << "# k=" << k << "\n";
    out << "# T_weight_units=" << horizon * unit << "\n";
    // The decoder's own `T_int`, read off the decoder rather than reconverted here — the §0 unit
    // rule. It is what the histogram bin width is a sixteenth of.
    out << "# T_time_units=" << horizon_time_units << "\n";
    out << "# mode=" << mode_name(mode) << "\n";
    out << "# shots=" << options.shots << "\n";
    out << "# seed=" << options.seed << "\n";
    out << "# bins_per_T=" << summary.weight_hist_bins_per_T << "\n";
    out << "# untimed=1 (the decomposition runs outside every timed window and is in no latency"
           " column)\n";
    out << "kind,key,value\n";

    auto scalar = [&](const char* key, double value) {
        out << "scalar," << key << "," << value << "\n";
    };
    scalar("shots_with_component_stats", (double)summary.shots_with_component_stats);
    scalar("shots_with_component_defects", (double)summary.shots_with_component_defects);
    scalar("mean_components", summary.mean_components);
    scalar("mean_singleton_components", summary.mean_singleton_components);
    scalar("mean_pair_components", summary.mean_pair_components);
    scalar("mean_nontrivial_components", summary.mean_nontrivial_components);
    scalar("mean_largest_component_size", summary.mean_largest_component_size);
    scalar("max_component_size", summary.max_component_size);
    scalar("mean_max_component_diameter_wint", summary.mean_max_component_diameter_wint);
    scalar("max_component_diameter_wint", summary.max_component_diameter_wint);
    scalar("boundary_touching_component_fraction", summary.boundary_touching_component_fraction);
    scalar("diameter_uncomputed_component_fraction", summary.diameter_uncomputed_component_fraction);
    scalar("frac_defects_in_trivial_components", summary.frac_defects_in_trivial_components);
    scalar("frac_defects_to_solver", summary.frac_defects_to_solver);
    scalar("frac_defects_committed_trivially", summary.frac_defects_committed_trivially);
    scalar("frac_defects_residual_trivially", summary.frac_defects_residual_trivially);
    scalar("solver_set_empty_rate", summary.solver_set_empty_rate);
    scalar("trivial_residual_rate", summary.trivial_residual_rate);
    // §A's run label. Over every shot, not only the analysed ones: the resolver ran on all of them.
    scalar("prune_component_max_size", (double)summary.prune_component_max_size);
    scalar("mean_defects_resolved_small", summary.mean_defects_resolved_small);

    auto histogram = [&](const char* kind, const std::vector<uint64_t>& bins) {
        for (size_t i = 0; i < bins.size(); i++)
            out << kind << "," << i << "," << bins[i] << "\n";
    };
    histogram("size_hist", summary.component_size_hist);
    histogram("diameter_hist", summary.component_diameter_hist);
    histogram("degree_hist", summary.h_degree_hist);
    histogram("edge_weight_hist", summary.h_edge_weight_hist);
    histogram("bcost_hist", summary.boundary_cost_hist);
}

/// One grid point: the per-shot rows, and the campaign accumulator the component file is written
/// from. `horizon_time_units` is the decoder's own `T_int`, carried so the component file can state
/// the unit its histogram bins are a sixteenth of without converting `T` a second time (§0).
struct PointResult {
    std::vector<ShotRow> rows;
    SpecMatchingAggregateStats stats;
    horizon_int horizon_time_units{0};
};

/// Both decoders on every shot, `H` first. Returns one row per shot, in shot order.
///
/// `progress` is ticked once per shot, warmup included, from between-shot positions only; see
/// `ProgressBar` for why that placement is the one that cannot disturb a measurement.
PointResult run_point(
    const stim::DetectorErrorModel& dem,
    const SpecMatchingConfig& config,
    const std::vector<std::vector<uint64_t>>& shots,
    size_t warmup,
    ProgressBar& progress) {
    auto decoder = SpecMatchingDecoder::from_detector_error_model(dem, config, NUM_DISTINCT_WEIGHTS);
    // Stock exact decode gets its own instance, so timing it does not disturb the front end's
    // `Mwpm(G)` — which the escalating path shares with the negative-weight preamble.
    pm::Mwpm reference = pm::detector_error_model_to_mwpm(dem, NUM_DISTINCT_WEIGHTS, false);
    std::vector<uint8_t> obs(std::max<size_t>(1, decoder.num_observables), 0);
    pm::total_weight_int weight = 0;

    size_t warmup_shots = std::min(warmup, shots.size());
    progress.begin(warmup_shots + shots.size());
    size_t retired = 0;
    for (size_t i = 0; i < warmup_shots; i++) {
        decoder.decode_to_obs(shots[i], obs.data(), weight, nullptr);
        pm::decode_detection_events(reference, shots[i], obs.data(), weight, false);
        progress.tick(++retired);
    }

    PointResult result;
    result.horizon_time_units = decoder.horizon;
    std::vector<ShotRow>& rows = result.rows;
    rows.reserve(shots.size());
    SpecMatchingProfile profile;
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

        // §C. Structural, and collected after this shot's `total_ns` was read, so nothing above is
        // charged for it. `decode_to_obs` does not accumulate — only `decode_batch` does — so the
        // campaign totals are gathered here, shot by shot, exactly as this driver gathers the rest.
        const ComponentStats& components = profile.components;
        row.h_defects = ball.h_nodes;
        row.components = components.num_components;
        row.singleton_components = components.num_singleton_components;
        row.pair_components = components.num_pair_components;
        row.nontrivial_components = components.num_nontrivial_components;
        row.largest_component = components.largest_component_size;
        row.defects_in_trivial = components.defects_in_trivial_components;
        row.defects_to_solver = components.defects_to_solver;
        // §A. Filled whether or not the component statistics were collected: it is the resolver's
        // own tally from the decode, not part of §C's post-shot decomposition.
        row.defects_resolved_small = profile.defects_resolved_small;
        result.stats.accumulate(profile);
        if (components.measured)
            result.stats.accumulate_component_histograms(decoder.component_histograms);

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
        // Both probe intervals are closed and both timers are stopped here, and the next shot's
        // decode opens its own probe inside `decode_to_obs`. A redraw that blocks on the terminal
        // therefore falls in the gap between two measurements rather than inside either.
        progress.tick(++retired);
    }
    progress.finish();
    return result;
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
    // The banner goes to stderr as well as being implied by the line above, because a wall-clock run
    // is a run whose every number is an upper bound and that fact must survive `> run.log`: stdout is
    // where the table goes and is routinely redirected, stderr is what stays on the screen.
    pm::perf::warn_if_wall_clock(stderr);
    const std::string probe = PreemptionProbe::mechanism_name();
    std::printf(
        "contamination probe: %s%s\n",
        probe.c_str(),
        probe == "none" ? " — no shot can be flagged on this platform; `contam` will read 0 because"
                          " nothing was measured, not because the machine was quiet"
                        : "");
    if (probe == "wall_minus_thread_cpu") {
        // The macOS path. Says what "contaminated" means here, since it is a threshold rather than an
        // exact counter, and where the threshold came from.
        std::printf(
            "  a shot is contaminated when wall time exceeds this thread's CPU time by more than"
            " %llu ns\n  (SPECMATCHING_PREEMPTION_SLACK_NS)\n",
            (unsigned long long)PreemptionProbe::slack_ns());
    }
    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }
    ProgressBar progress;
    progress.interactive = stderr_is_terminal();
    progress.enabled = options.progress == ProgressMode::ALWAYS ||
                       (options.progress == ProgressMode::AUTO && progress.interactive);
    size_t jobs = options.distances.size() * options.error_rates.size() * options.horizons.size() *
                  options.modes.size() * options.ks.size();
    size_t job = 0;

    std::printf(
        "writing one latency log%s per (d, p, T, mode, k) to %s\n\n",
        options.component_stats ? " and one component file" : "",
        options.out_dir.c_str());
    // `sparse_mean` is the serial machine (H, then G on a rejected shot) and `system_mean` the
    // concurrent one (H and G at once, ending at the first usable matching); `speedup` is
    // `stock_mean / system_mean` and nothing else. `b+h` is §D's figure of merit — the solve on `H`
    // plus the harvest, the two stages a smaller node set moves — `triv`/`solver` are the share of
    // the defect set that sits in a component of size <= 2 and the share the solver would still have
    // to take, and `resolved` is the mean number of defects the resolver took off the solver per shot
    // at this run's `k`.
    std::printf(
        "%4s %8s %5s %8s %3s %8s %11s %11s %11s %8s %11s %8s %9s %8s %6s %6s %9s\n",
        "d",
        "p",
        "T",
        "mode",
        "k",
        "shots",
        "stock_mean",
        "sparse_mean",
        "system_mean",
        "speedup",
        "b+h_mean",
        "escal",
        "escal_frac",
        "contam",
        "triv",
        "solver",
        "resolved");

    size_t files_written = 0;
    for (size_t distance : options.distances) {
        for (double noise : options.error_rates) {
            {
                char sampling[128];
                std::snprintf(
                    sampling, sizeof(sampling), "sampling %zu shots at d=%zu p=%g ...", options.shots, distance, noise);
                progress.note(sampling);
            }
            Experiment experiment = generate(distance, distance, noise, options.shots, options.seed);
            pm::Mwpm probe = pm::detector_error_model_to_mwpm(experiment.dem, NUM_DISTINCT_WEIGHTS, false);
            double unit = edge_weight_units(probe.flooder.graph);

            for (double horizon : options.horizons) {
                for (BallGraphBuildMode mode : options.modes) {
                    // One independent experiment per `k`, in the order the flag listed them. The
                    // corpus, the seed and the warmup are the same for all of them — only `k`
                    // differs — so the runs are read against each other directly, and `k = 0` is
                    // the un-pruned side when it is in the list.
                    for (int k : options.ks) {
                        char label[128];
                        std::snprintf(
                            label,
                            sizeof(label),
                            "[%zu/%zu] d=%zu p=%g T=%g %s k=%d",
                            ++job,
                            jobs,
                            distance,
                            noise,
                            horizon,
                            mode_name(mode),
                            k);
                        progress.label = label;
                        PointResult point = run_point(
                            experiment.dem,
                            config_for(horizon, unit, mode, options.component_stats, k),
                            experiment.shots,
                            options.warmup,
                            progress);
                        const std::vector<ShotRow>& rows = point.rows;

                        std::string path =
                            options.out_dir + "/" +
                            point_name("latency", distance, noise, horizon, mode, k, options.tag);
                        std::ofstream out(path);
                        if (!out.is_open()) {
                            std::cerr << "error: could not open " << path << " for writing"
                                      << " (does " << options.out_dir << " exist?)\n";
                            return 1;
                        }
                        write_header(out, options, distance, noise, horizon, mode, k, unit, progress);
                        for (size_t i = 0; i < rows.size(); i++) {
                            const ShotRow& row = rows[i];
                            out << i << "," << row.defects << "," << row.escalated << "," << row.certified << ","
                                << row.contaminated << "," << row.stock_g_ns << "," << row.blossom_ns << ","
                                << row.dscan_ns << "," << row.hrvst_ns << "," << row.escal_stock_ns << ","
                                << row.excluded_ns << "," << row.total_ns << "," << row.h_defects << ","
                                << row.components << ","
                                << row.singleton_components << "," << row.pair_components << ","
                                << row.nontrivial_components << "," << row.largest_component << ","
                                << row.defects_in_trivial << "," << row.defects_to_solver << ","
                                << row.defects_resolved_small << "\n";
                        }
                        out.close();
                        files_written++;

                        // The component structure goes in its own file: it is per component and per
                        // edge rather than per shot, so it does not fit the latency schema, and
                        // keeping the two apart is what stops a distribution over components being
                        // read as a distribution over shots.
                        SpecMatchingSummary component_summary = summarize(point.stats);
                        if (options.component_stats) {
                            std::string component_path =
                                options.out_dir + "/" +
                                point_name("components", distance, noise, horizon, mode, k, options.tag);
                            std::ofstream component_out(component_path);
                            if (!component_out.is_open()) {
                                std::cerr << "error: could not open " << component_path << " for writing\n";
                                return 1;
                            }
                            write_component_file(
                                component_out,
                                options,
                                distance,
                                noise,
                                horizon,
                                mode,
                                k,
                                unit,
                                point.horizon_time_units,
                                component_summary);
                            component_out.close();
                            files_written++;
                        }

                        QuickStats stock = quick_stats(rows, [](const ShotRow& row) { return row.stock_g_ns; });
                        QuickStats sparse = quick_stats(rows, [](const ShotRow& row) { return row.sparse_ns(); });
                        QuickStats system = quick_stats(rows, [](const ShotRow& row) { return row.system_ns(); });
                        QuickStats solve_harvest =
                            quick_stats(rows, [](const ShotRow& row) { return row.solve_and_harvest_ns(); });
                        size_t escalated = 0;
                        size_t contaminated = 0;
                        for (const ShotRow& row : rows) {
                            escalated += (size_t)row.escalated;
                            contaminated += (size_t)row.contaminated;
                        }
                        // The decoder system without sparsification over the same system predicating
                        // on `H`, both including every escalated shot: what the front end buys once
                        // the escalations are priced as the lost race they are.
                        double speedup = system.mean > 0 ? stock.mean / system.mean : 0.0;
                        double escalated_fraction = rows.empty() ? 0.0 : (double)escalated / (double)rows.size();
                        std::printf(
                            "%4zu %8g %5g %8s %3d %8zu %9.3fus %9.3fus %9.3fus %8.3fx %9.3fus %8zu %8.3f%% %8zu"
                            " %5.1f%% %5.1f%% %9.3f\n",
                            distance,
                            noise,
                            horizon,
                            mode_name(mode),
                            k,
                            system.kept,
                            stock.mean / 1000.0,
                            sparse.mean / 1000.0,
                            system.mean / 1000.0,
                            speedup,
                            solve_harvest.mean / 1000.0,
                            escalated,
                            100.0 * escalated_fraction,
                            contaminated,
                            100.0 * component_summary.frac_defects_in_trivial_components,
                            100.0 * component_summary.frac_defects_to_solver,
                            component_summary.mean_defects_resolved_small);
                        std::fflush(stdout);
                    }
                }
            }
        }
    }

    std::printf(
        "\n%zu file(s) written. Means above are over uncontaminated shots only, escalated shots"
        " included;\n`contam` is how many of the %zu shots at each point were dropped from them, and"
        " `escal`/`escal_frac`\nare over all %zu.\n"
        "`sparse_mean` is the serial machine: blossom+dscan+hrvst, plus the Phase-2 re-decode after it"
        " on the\nshots that escalate. Ball intersect, H build and Mwpm(H) build are in neither"
        " column and are logged\nas `excluded_ns` in every file.\n"
        "`system_mean` is the concurrent decoder system — H and G solved at once, the shot ending at"
        " the first\nusable matching. An escalating shot is charged the stock decode on G alone,"
        " because that decode was\nalready running and H's work bought nothing; any other shot is"
        " charged min(blossom+dscan+hrvst,\nstock). Escalation is priced as a lost race rather than"
        " as an added tail.\n"
        "`speedup` is stock_mean/system_mean: the decoder system without graph sparsification against"
        " the same\nsystem predicating on H, over the same shots. `sparse_mean` is never a"
        " denominator.\n"
        "`b+h_mean` is blossom+hrvst alone — the two stages a smaller node set moves. `triv` and"
        " `solver`\nare shares of H's defects: in a component of size <= 2, and left to the solver"
        " anyway. Both are\nstructural, measured outside every timed window, and detailed per grid"
        " point in components_*.csv.\n"
        "Each row is one k: the solver saw only components of size > k, and `resolved` is the mean"
        " number of\ndefects the resolver settled off it per shot. That resolve runs in series"
        " before the solve, on the\ncritical path, and is in none of the times above — this driver"
        " measures solver + harvest on the\nsize-> k graph, and the resolve is outside that scope by"
        " intent, to be measured separately.\n",
        files_written,
        options.shots,
        options.shots);
    return 0;
}
