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

/// The per-component run-time profiler.
///
/// For the `d x d x d` rotated surface-code memory experiment under circuit-level noise at rate
/// `p`, at horizons `T`, this logs the run time of truncated sparse blossom on **every connected
/// component of `H`**, keyed by component size. One row per component. It performs no analysis of
/// that table and prints no interpretation of it.
///
/// Companion to `sparse_graph_stats`, which measures the same components' *structure* and does so
/// with no clock at all. This one measures nothing but their cost.
///
/// ## What is timed (§1)
///
/// One timer pair per component, around the existing per-component path and nothing else:
///
/// ```
///   t0 = hires_now_ns()
///     build_component_subgraph(H, split, c)          // sub-H(C)              §2.4
///     h_mwpm.rebuild(tables, sub, arena)             // the instance on it    §2.3
///     process_timeline_until_horizon(instance, dets_C, T_int)
///     if COMPLETE: extract_only_to_obs(...)          // extract obs and weight
///   t1 = hires_now_ns()
///   if TRUNCATED: abandon_shot(instance)             // instance.reset(), NOT timed
/// ```
///
/// Recorded with `t1 - t0`: `comp_size = |C|` and `status`. **Status is kept** because a truncated
/// component runs to the horizon and its cost at a given size is a different quantity from a
/// completed one's; the two cannot be separated after the fact.
///
/// Components of size `<= 2` get **no row**. They are still built, solved and abandoned on the same
/// path — the workload is untouched — but a solve that small is of the same order as the clock reads
/// bracketing it, so the row would report the timer rather than the component. `run.log` records the
/// cut as `min_logged_comp_size`; `components.csv` therefore holds no `comp_size` below it, and a
/// consumer counting components per shot must not read it off this table.
///
/// The stock escalation is not run and not timed — `escalate_to_stock` is never called and this
/// file does not include it. The profiler stops at the trigger. Observables and weights are
/// discarded on every component.
///
/// ## Why the loop is here and not a flag on the decoder (§6.1, §6.2)
///
/// Hard constraint 1 is **no decoder change**, and hard constraint 2 is **one timer pair per
/// component**. Those two together decide the shape of this file. `BallDecoder::decode_impl` runs
/// the per-component loop, but the only way to get a number out of it is to pass a `BallProfile`,
/// and a profiled shot arms the decoder's own three stage timers *inside* the region §1 asks for a
/// single pair around — six extra readings per component, on a region whose smallest rows are the
/// same order as one reading. It also sums its stages over the shot rather than reporting them per
/// component, which is the quantity this binary exists to produce.
///
/// So the loop is driven from here, out of the same primitives, with `prof == nullptr` throughout
/// and exactly two clock reads per component. It is `decode_impl`'s loop restricted to the
/// production path — the same `build_ball_graph`, the same `decompose_ball_components`, the same
/// `build_component_subgraph`, the same `BallMwpm::rebuild`, the same
/// `process_timeline_until_horizon`, the same `extract_only_to_obs`, the same `abandon_shot` — and
/// it drops only what §1 puts outside the region and what a discarded output does not need:
/// `accumulate_component_harvest`, the residual remap (empty on the production path), the
/// `sort_residual` pass, the profile fills and the `verify_*` hooks. Nothing under
/// `src/specmatching/` is touched.
///
/// `BallDecoder::truncated_scheme_escalates` is the same pattern inside the decoder itself: the
/// per-component primitives driven directly by a caller that wants one specific thing off them.
///
/// ## Units
///
/// `--T` is in **multiples of one lattice edge weight**, exactly as in `sparse_graph_stats` — the
/// normalisation under which the `M1` structural results are `d`-independent, and the one every
/// horizon in `benchmarks/spec_matching/` is quoted in. That is what makes `T = 2` here the same
/// `T = 2` as in the companion's tables, which is the whole point of the two binaries sharing a
/// grid. The design's §3 quotes `--T` in raw DEM float weight units; that reading is a deviation
/// recorded rather than silently taken, and it is the same one `sparse_graph_stats` records, for
/// the same reason: a horizon of 2 raw DEM float units is a fraction of an edge at `p = 1e-3` and
/// would truncate every component. All three readings — `T` (edge weights), `T_weight_units` (DEM
/// float) and `T_int` (the flooder's integer time units, read off the decoder rather than
/// reconverted) — go into `run.log`, so no consumer has to infer which is which.
///
/// The file also lives in `benchmarks/spec_matching/` rather than the design's
/// `benchmarks/two_phase/`, which is this tree's name for that corpus; the generator call is the
/// corpus's own (see `generator_call` below for its two deviations, both inherited).
///
/// ## Usage
///
///   component_profiler [--d 5,7,9] [--p 1e-3,5e-4] [--T 1.5,2] [--shots N] [--warmup 1000]
///                      [--seed S] [--out DIR]
///
/// Writes `components.csv` and `run.log` to `--out`. The agent runs no campaigns (§6.7); smoke
/// runs at `--shots 100 --d 5` exist only to check that the output parses.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmarks/spec_matching/hires_timer.h"
#include "benchmarks/spec_matching/profiler_util.h"
#include "specmatching/spec_matching/manifold/ball_decoding.h"

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#ifndef SPECMATCHING_GIT_HASH
#define SPECMATCHING_GIT_HASH "unknown"
#endif
#ifndef SPECMATCHING_BUILD_FLAGS
#define SPECMATCHING_BUILD_FLAGS "unknown"
#endif
#ifndef SPECMATCHING_STIM_VERSION
#define SPECMATCHING_STIM_VERSION "unknown"
#endif

using namespace pm::spec_matching;
using namespace pm::spec_matching::profiler;

namespace {

/// Shots sampled per batch. Bounds the memory a large cell needs — a `d = 17` shot's detection
/// events cost most of a kilobyte — without materialising the corpus. The stream is a function of
/// the seed and the batch sizes alone, and both are identical across the `T` sweep, which is what
/// makes shot `i` of `(d, p)` the same shot at every horizon (§3's `--seed`).
constexpr size_t SHOT_BATCH = 4096;

/// §5's fixed row arena. Reserved once at startup and never grown, so appending a row allocates
/// nothing (§6.6); drained to `components.csv` between batches and at the end of every cell, from
/// between-shot positions only — never from inside a timed region.
///
/// A shot cannot overflow it: a component needs a defect, so a shot has at most `|H|` components,
/// and `|H|` is bounded by the detector count of the largest circuit in the grid — tens of
/// thousands, three orders below the flush threshold.
constexpr size_t ROW_ARENA_CAPACITY = 1u << 20;
constexpr size_t ROW_ARENA_FLUSH_AT = ROW_ARENA_CAPACITY / 2;

/// The smallest component that gets a row. A size-1 or size-2 solve is of the same order as the two
/// clock reads around it, so its row carries more timer than component; those sizes are dropped at
/// the logging point rather than kept and filtered downstream.
///
/// Only the row is dropped. Every component of every shot is still built, solved to the horizon and
/// abandoned exactly as before — the workload a measured row sits in the middle of is unchanged, and
/// so is the state each component leaves for the next one.
constexpr uint32_t MIN_LOGGED_COMP_SIZE = 3;

struct Options {
    std::vector<size_t> distances = {5};
    std::vector<double> error_rates = {1e-3};
    std::vector<double> horizons = {1.5, 2.0};
    size_t shots = 10000;
    size_t warmup = 1000;
    uint64_t seed = 20260905;
    std::string out_dir = "benchmarks/spec_matching/results/component_profiler";
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
        if (flag == "--d" || flag == "--distances") {
            options.distances = parse_list<size_t>(next());
        } else if (flag == "--p" || flag == "--error-rates") {
            options.error_rates = parse_list<double>(next());
        } else if (flag == "--T" || flag == "--horizons") {
            options.horizons = parse_list<double>(next());
        } else if (flag == "--shots") {
            options.shots = std::stoul(next());
        } else if (flag == "--warmup") {
            options.warmup = std::stoul(next());
        } else if (flag == "--seed") {
            options.seed = std::stoull(next());
        } else if (flag == "--out" || flag == "--out-dir") {
            options.out_dir = next();
        } else {
            throw std::invalid_argument("unrecognised flag " + flag);
        }
    }
    if (options.distances.empty() || options.error_rates.empty() || options.horizons.empty())
        throw std::invalid_argument("--d, --p and --T each need at least one value");
    if (options.shots == 0)
        throw std::invalid_argument("--shots must be positive");
    for (size_t distance : options.distances) {
        if (distance < 3 || distance % 2 == 0)
            throw std::invalid_argument("--d values must be odd and at least 3");
    }
    for (double noise : options.error_rates) {
        if (!(noise > 0) || !(noise < 1))
            throw std::invalid_argument("--p values must lie in (0, 1)");
    }
    for (double horizon : options.horizons) {
        if (!(horizon > 0))
            throw std::invalid_argument("--T values must be positive");
    }
    return options;
}

/// One row of `components.csv`: one component of one shot.
///
/// `d`, `p` and `T` are cell constants and are written at flush time rather than stored per row —
/// at the largest `(d, p)` a shot has hundreds of components, so the arena is the one structure
/// here whose per-row width is worth thinking about.
struct Row {
    uint64_t shot;
    uint32_t comp_size;
    /// `1` iff the component's own truncated solve did not finish by `T`.
    uint8_t truncated;
    /// The raw, **uncorrected** `hires_now_ns()` delta. The analysis subtracts `timer_overhead_ns`.
    uint64_t ns;
};

/// The cell labels a flush stamps onto every row it writes.
struct CellLabels {
    size_t distance{0};
    double noise{0};
    double horizon{0};
};

struct RowArena {
    std::vector<Row> rows;
    CellLabels labels;
    std::ofstream* out{nullptr};
    /// Rows written this cell, and over the whole run.
    uint64_t rows_this_cell{0};
    uint64_t rows_total{0};

    void configure(std::ofstream& stream) {
        out = &stream;
        rows.reserve(ROW_ARENA_CAPACITY);
    }

    inline void append(uint64_t shot, uint32_t comp_size, bool truncated, uint64_t ns) {
        rows.push_back(Row{shot, comp_size, (uint8_t)(truncated ? 1 : 0), ns});
    }

    /// Called from between-shot positions only. Never from inside a timed region: the whole point
    /// of the arena is that the write to the file happens somewhere no measurement can see it.
    void flush() {
        for (const Row& row : rows) {
            *out << labels.distance << "," << labels.noise << "," << labels.horizon << "," << row.shot << ","
                 << row.comp_size << "," << (row.truncated ? "TRUNCATED" : "COMPLETE") << "," << row.ns << "\n";
        }
        rows_this_cell += rows.size();
        rows_total += rows.size();
        // `clear` on a vector of trivial rows keeps the capacity, so the arena is reserved once and
        // reused for the life of the process.
        rows.clear();
        out->flush();
    }
};

/// Whether the progress line may use carriage returns and the erase-to-end-of-line escape. A
/// redirected run gets one plain line per batch instead of a megabyte of control codes.
bool stderr_is_terminal() {
#if defined(__unix__) || defined(__APPLE__)
    if (isatty(fileno(stderr)) != 1)
        return false;
    const char* term = std::getenv("TERM");
    return term != nullptr && std::string(term) != "dumb";
#else
    return false;
#endif
}

constexpr size_t PROGRESS_BAR_WIDTH = 24;

/// The campaign's progress: a bar over the current cell's shots, and the cell's place in the grid.
///
/// **No ETA and no shots-per-second**, for the same reason `sparse_graph_stats` gives: every number
/// on the line is a ratio of counts the loop already keeps, so drawing it reads no clock. This
/// binary does own a clock, but it owns it for the rows, and a progress line is not a measurement.
/// Drawn between batches, to stderr, from outside every timed region.
void draw_progress(
    bool interactive, const char* label, size_t cell_index, size_t total_cells, size_t done, size_t cell_shots) {
    double fraction = cell_shots ? (double)done / (double)cell_shots : 1.0;
    if (!interactive) {
        std::fprintf(
            stderr,
            "[%zu/%zu] %s  %5.1f%%  %zu/%zu shots\n",
            cell_index,
            total_cells,
            label,
            100.0 * fraction,
            done,
            cell_shots);
        std::fflush(stderr);
        return;
    }
    char bar[PROGRESS_BAR_WIDTH + 1];
    size_t filled = (size_t)(fraction * (double)PROGRESS_BAR_WIDTH);
    for (size_t i = 0; i < PROGRESS_BAR_WIDTH; i++)
        bar[i] = i < filled ? '#' : '.';
    bar[PROGRESS_BAR_WIDTH] = '\0';
    std::fprintf(
        stderr,
        "\r[%zu/%zu] %s  [%s] %5.1f%%  %zu/%zu shots\x1b[K",
        cell_index,
        total_cells,
        label,
        bar,
        100.0 * fraction,
        done,
        cell_shots);
    std::fflush(stderr);
}

/// What stands in for the bar while a `(d, p)` is being set up: the DEM and the ball tables are
/// built before the first shot is sampled, and at `d = 17` and up that is long enough to look like
/// a hang.
void draw_setup(bool interactive, const char* label, size_t cell_index, size_t total_cells) {
    std::fprintf(
        stderr,
        interactive ? "\r[%zu/%zu] %s  building DEM and ball tables...\x1b[K"
                    : "[%zu/%zu] %s  building DEM and ball tables...\n",
        cell_index,
        total_cells,
        label);
    std::fflush(stderr);
}

/// The generator call, as one string, so `run.log` records what `(d, p)` actually meant.
///
/// It is `profiler_util.h`'s `ShotSampler::make` — the same call the rest of the
/// `benchmarks/spec_matching/` corpus is generated by, which is the requirement: `(d, p)` has to
/// mean here what it means in every table beside it, `sparse_graph_stats`'s included.
///
/// Two deviations from the design's prose, recorded rather than silently taken, and both inherited
/// from that shared call: the basis is `rotated_memory_x`, not `_z`; and three noise parameters are
/// set, not four (`before_round_data_depolarization` is left at 0).
std::string generator_call(size_t distance, size_t rounds, double noise) {
    std::ostringstream out;
    out << "stim::CircuitGenParameters(rounds=" << rounds << ", distance=" << distance
        << ", task=\"rotated_memory_x\"); after_clifford_depolarization=" << noise
        << "; after_reset_flip_probability=" << noise << "; before_measure_flip_probability=" << noise
        << "; before_round_data_depolarization=0 (unset)"
        << "; stim::generate_surface_code_circuit(...).circuit"
        << "; stim::ErrorAnalyzer::circuit_to_detector_error_model(circuit, false, true, false, 0, false, false)";
    return out.str();
}

/// The three "did a buffer have to grow?" counters, which is how invariant 11 is stated in this
/// tree. Read at the end of the warm-up and again at the end of the cell; if they moved, a measured
/// shot allocated.
struct GrowEvents {
    uint64_t arena{0};
    uint64_t split{0};
    uint64_t mwpm{0};

    static GrowEvents of(const BallDecoder& decoder) {
        return GrowEvents{decoder.arena.grow_events, decoder.arena.split.grow_events, decoder.h_mwpm.grow_events};
    }
    bool operator==(const GrowEvents& rhs) const {
        return arena == rhs.arena && split == rhs.split && mwpm == rhs.mwpm;
    }
};

/// §1's measurement, and the whole of it.
///
/// Everything before the loop is `decode_impl`'s untimed preamble: the negative-weight preamble on
/// `G`, the ball intersection that builds `H`, and the union-find that splits it. None of it is in
/// any row — §1's region starts at sub-`H` construction for a component.
///
/// `record` says whether the shot's rows are kept; a warm-up shot runs the identical path and
/// throws them away, so it costs the same and touches the same buffers as a measured one.
///
/// `part` is the caller's reused `HarvestResult`, assigned into exactly as `decode_impl` assigns
/// into its own, so the two paths have the same allocation behaviour. Its contents are never read:
/// this binary keeps no decoder output.
void profile_shot(
    BallDecoder& decoder,
    const std::vector<uint64_t>& shot,
    uint64_t shot_index,
    bool record,
    HarvestResult& part,
    RowArena& arena) {
    decoder.compute_seeded_detection_events(shot, decoder.seeded_scratch);
    build_ball_graph(
        decoder.tables, decoder.seeded_scratch, decoder.horizon, decoder.arena, decoder.config.mode, nullptr, nullptr);
    const BallGraph& h = decoder.arena.graph;
    BallComponentSplit& split = decoder.arena.split;
    decompose_ball_components(h, split);

    // A zero-defect shot decomposes into no components and therefore produces no rows (§4).
    for (size_t c = 0; c < split.num_components(); c++) {
        // ---- §1. Two reads, and nothing between them but the per-component path.
        uint64_t started = hires_now_ns();
        build_component_subgraph(h, split, c);
        decoder.h_mwpm.rebuild(decoder.tables, split.sub, decoder.arena, nullptr);
        // The component's nodes are `0..s-1` by construction, so its detection events are every
        // node. `h_dets_scratch` keeps its capacity across components and shots.
        decoder.h_dets_scratch.clear();
        decoder.h_dets_scratch.reserve(split.sub.num_nodes());
        for (size_t i = 0; i < split.sub.num_nodes(); i++)
            decoder.h_dets_scratch.push_back(i);
        TimelineStatus status =
            process_timeline_until_horizon(decoder.h_mwpm.mwpm, decoder.h_dets_scratch, decoder.horizon);
        if (status == TimelineStatus::COMPLETE) {
            // §M3.4's production extraction: the observables and the weight, no harvest.
            part = decoder.harvester.extract_only_to_obs(decoder.h_mwpm.mwpm, decoder.h_dets_scratch);
        }
        uint64_t ended = hires_now_ns();

        // ---- Outside the region. `abandon_shot` is `Mwpm::reset`, which §1 puts outside the timed
        // window; a `COMPLETE` component was already left clean by the extraction's own
        // `reset_for_next_shot`.
        if (status == TimelineStatus::TRUNCATED)
            abandon_shot(decoder.h_mwpm.mwpm);

        split.status[c] = status == TimelineStatus::TRUNCATED ? BallComponentSplit::COMPONENT_TRUNCATED
                                                              : BallComponentSplit::COMPONENT_COMPLETE;
        if (record && split.sizes[c] >= MIN_LOGGED_COMP_SIZE)
            arena.append(shot_index, split.sizes[c], status == TimelineStatus::TRUNCATED, ended - started);
    }
}

std::string stamp() {
    std::time_t now = std::time(nullptr);
    char buffer[64];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", std::localtime(&now));
    return std::string(buffer);
}

/// One `(d, p)` of the grid: circuit, DEM, ball tables, then the `T` sweep over the same shots.
///
/// Throws on anything that makes the cell unmeasurable; `main` logs it and stops. `cell_index` is
/// the caller's running position in the `(d, p, T)` grid, advanced here.
void run_dp_cell(
    const Options& options,
    size_t distance,
    double noise,
    RowArena& arena,
    std::ofstream& log,
    bool interactive,
    size_t& cell_index,
    size_t total_cells) {
    char setup_label[160];
    std::snprintf(setup_label, sizeof(setup_label), "d=%zu p=%g", distance, noise);
    draw_setup(interactive, setup_label, cell_index + 1, total_cells);

    // ---- §4 step 1. Circuit -> DEM, through the corpus's own generator call. The negative-weight
    // preamble is the decoder's, and `dem_has_negative_weights` below records whether it had
    // anything to do (§M2.1).
    ShotSampler sampler = ShotSampler::make(distance, distance, noise, options.seed);
    log << "\ncorpus d=" << distance << " rounds=" << distance << " p=" << noise << "\n";
    log << "  generator_call=" << generator_call(distance, distance, noise) << "\n";

    // ---- §4 step 2. Ball tables once, at `T_max = max(T list)` and `R = 2 * T_max`.
    pm::Mwpm probe = pm::detector_error_model_to_mwpm(sampler.dem, NUM_DISTINCT_WEIGHTS, false);
    double unit = edge_weight_units(probe.flooder.graph);
    double max_horizon = *std::max_element(options.horizons.begin(), options.horizons.end());

    BallConfig ball_config;
    ball_config.T = max_horizon * unit;
    ball_config.ball.T_max = max_horizon * unit;
    ball_config.ball.R = 2.0 * ball_config.ball.T_max;
    // Nothing in this binary reads the component structure — that is the companion's job — and
    // nothing in it verifies, escalates, or collects a stage split.
    ball_config.collect_component_stats = false;
    ball_config.collect_harvest_diagnostics = false;
    ball_config.collect_structural_counters = false;
    ball_config.verify_component_decomposition = false;
    ball_config.verify_against_g = false;

    // §4 step 2's assert. `R >= 2 * T_max` is what makes a pair within `2T` reachable at all, and
    // an undersized `R` changes the answer without any crash, so it is checked in release builds
    // too — `validate()` is the decoder's own version of the same check.
    ball_config.ball.validate();
    if (!(ball_config.ball.R >= 2.0 * ball_config.ball.T_max))
        throw std::runtime_error("the compiled ball radius is below 2 * T_max");

    BallDecoder decoder = BallDecoder::from_detector_error_model(sampler.dem, ball_config, NUM_DISTINCT_WEIGHTS);
    // Set explicitly rather than left at their defaults: `decode_impl` writes both from the config
    // on every shot, and this loop is standing in for it.
    decoder.harvester.collect_diagnostics = false;
    decoder.harvester.use_legacy_enumeration = false;
    log << "  dem_has_negative_weights=" << (decoder.dem_has_negative_weights ? 1 : 0)
        << " edge_weight_unit=" << unit << " T_max=" << ball_config.ball.T_max << " R=" << ball_config.ball.R << "\n";
    log.flush();

    HarvestResult part;

    // ---- §4 steps 3-4. One `T` at a time, over the same shots: the sampler's stream is a function
    // of its rng and of the batch sizes, and both are reset identically below, so shot `i` of
    // `(d, p)` is the same shot at every horizon. The tables are not recompiled — `set_horizon` is
    // exactly the retarget §4 step 2 leaves room for.
    for (double horizon : options.horizons) {
        cell_index++;
        char label[160];
        std::snprintf(label, sizeof(label), "d=%zu p=%g T=%g", distance, noise, horizon);

        decoder.set_horizon(horizon * unit);
        sampler.rng = std::mt19937_64(options.seed);

        arena.labels = CellLabels{distance, noise, horizon};
        arena.rows_this_cell = 0;
        log << "  cell T=" << horizon << " T_weight_units=" << horizon * unit
            << " T_int=" << (int64_t)decoder.horizon << " started=" << stamp() << "\n";
        log.flush();

        size_t done = 0;
        size_t cell_total = options.warmup + options.shots;
        // Read before the first shot so that `--warmup 0` compares against the decoder's
        // as-constructed state rather than against zeros.
        GrowEvents after_warmup = GrowEvents::of(decoder);
        draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
        while (done < cell_total) {
            size_t batch = std::min(SHOT_BATCH, cell_total - done);
            sampler.sample(batch);
            for (size_t i = 0; i < batch; i++) {
                size_t index = done + i;
                bool measured = index >= options.warmup;
                // `shot` is the index within the cell, 0-based **after** warm-up.
                uint64_t shot_index = measured ? (uint64_t)(index - options.warmup) : 0;
                profile_shot(decoder, sampler.shots[i], shot_index, measured, part, arena);
                if (index + 1 == options.warmup)
                    after_warmup = GrowEvents::of(decoder);
            }
            done += batch;
            // Between batches, and from a between-shot position: no timed region is open, so
            // neither the flush nor the progress line lands in any row.
            if (arena.rows.size() >= ROW_ARENA_FLUSH_AT)
                arena.flush();
            draw_progress(interactive, label, cell_index, total_cells, done, cell_total);
        }
        arena.flush();
        if (interactive)
            std::fprintf(stderr, "\r\x1b[K");

        // §6.6. Invariant 11 as this tree states it: after the warm-up no buffer grew, so no
        // measured shot allocated. Reported rather than asserted, because a cell whose warm-up was
        // too short for its corpus should say so in its own artifact.
        GrowEvents at_end = GrowEvents::of(decoder);
        log << "  cell T=" << horizon << " finished=" << stamp() << " rows=" << arena.rows_this_cell
            << " no_growth_after_warmup=" << ((at_end == after_warmup) ? 1 : 0) << " (grow_events arena/split/mwpm "
            << after_warmup.arena << "/" << after_warmup.split << "/" << after_warmup.mwpm << " -> " << at_end.arena
            << "/" << at_end.split << "/" << at_end.mwpm << ")\n";
        log.flush();
    }
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

    std::error_code dir_error;
    std::filesystem::create_directories(options.out_dir, dir_error);
    if (dir_error) {
        std::cerr << "error: could not create " << options.out_dir << ": " << dir_error.message() << "\n";
        return 1;
    }

    std::ofstream components(options.out_dir + "/components.csv");
    std::ofstream log(options.out_dir + "/run.log");
    if (!components.is_open() || !log.is_open()) {
        std::cerr << "error: could not open the output files under " << options.out_dir << "\n";
        return 1;
    }
    components << "d,p,T,shot,comp_size,status,ns\n";

    RowArena arena;
    arena.configure(components);

    // ---- §2's startup calibration, before any decoding, so it measures the clock rather than the
    // clock plus whatever the decoder left in the caches.
    double overhead_ns = hires_timer_overhead_ns();
    bool thread_scoped = hires_timer_is_thread_scoped();

    log << "component_profiler - per-component run time of truncated sparse blossom on H\n";
    log << "git_hash=" << SPECMATCHING_GIT_HASH << "\n";
    log << "build_flags=" << SPECMATCHING_BUILD_FLAGS << "\n";
    log << "stim_version=" << SPECMATCHING_STIM_VERSION << "\n";
#if defined(__APPLE__)
    log << "platform=macos\n";
#else
    log << "platform=linux\n";
#endif
    log << "hires_timer_name=" << hires_timer_name() << "\n";
    log << "timer_is_thread_scoped=" << (thread_scoped ? 1 : 0) << "\n";
    log << "timer_overhead_ns=" << overhead_ns << " (10^6 back-to-back hires_now_ns() calls)\n";
    log << "rows_are_uncorrected=1 (ns is the raw delta; the analysis subtracts timer_overhead_ns)\n";
    log << "min_logged_comp_size=" << MIN_LOGGED_COMP_SIZE
        << " (components of size <= 2 are solved and timed as before but get no row: the timer overhead is of"
           " their own order, so their times are overhead rather than measurement)\n";
    log << "seed=" << options.seed << "\n";
    log << "shots_per_cell=" << options.shots << " (after warm-up)\n";
    log << "warmup_shots_per_cell=" << options.warmup << " (run and discarded)\n";
    log << "threads=1\n";
    log << "T_unit=one lattice edge weight (the median discretised edge of G, in DEM float units);"
           " each cell below records T, T_weight_units and T_int\n";
    log << "escalation=not run and not timed; escalate_to_stock is never called\n";

    // A wall-clock fallback does not fail and does not look wrong — it just quietly charges every
    // component for whatever else the machine was doing — so it is said loudly in both places.
    if (!thread_scoped) {
        log << "WARNING: no thread-scoped clock was available, so every ns on this run includes time"
               " the thread spent off the CPU and is an upper bound only.\n";
        std::fprintf(
            stderr,
            "\n  WARNING: component times are being measured with a WALL CLOCK, backend '%s'.\n"
            "  Every ns includes time the thread spent off the CPU and is an upper bound only.\n\n",
            hires_timer_name());
        std::fflush(stderr);
    }
    log.flush();

    bool interactive = stderr_is_terminal();
    size_t total_cells = options.distances.size() * options.error_rates.size() * options.horizons.size();
    size_t cell_index = 0;
    bool ok = true;

    for (size_t distance : options.distances) {
        if (!ok)
            break;
        for (double noise : options.error_rates) {
            try {
                run_dp_cell(options, distance, noise, arena, log, interactive, cell_index, total_cells);
            } catch (const std::exception& error) {
                std::fprintf(stderr, "\n");
                std::cerr << "error at d=" << distance << " p=" << noise << ": " << error.what() << "\n";
                log << "  ERROR=" << error.what() << "\n";
                ok = false;
                break;
            }
        }
    }

    log << "\nfinished=" << stamp() << "\n";
    log << "rows_total=" << arena.rows_total << "\n";
    log << (ok ? "status=ok\n" : "status=incomplete (see ERROR above)\n");

    std::printf("wrote %s/{components.csv,run.log}\n", options.out_dir.c_str());
    return ok ? 0 : 1;
}
