# Latency speedups

Per-shot decode latency of the spec-matching decoder against stock exact sparse blossom, over a
72-point sweep of a rotated surface-code memory-X experiment, one million shots a point. This file
is the record of that measurement.

It replaces an earlier report written from `profiler_driver`, which no longer exists. That binary
measured a monolithic solve on `H` behind a lookup-table small-component resolver, and both are
gone: every component of `H` is now decided by truncated sparse blossom run on it in isolation, and
the components are spread over `k` solver cores by a load balancer. Nothing from the old tables
carries over — different decomposition, different input to the solver, different definition of the
timed region — so none of it is reproduced here.

The numbers below are a **latency** read of a **modelled** `k`-core critical path. Read the model
and its assumptions before the tables; the assumptions are load-bearing, not boilerplate.

## What was measured

Per shot the profiler logs three latencies. `run.log` in the run directory carries the definitions
verbatim; they are:

* **`fallback`** — stock exact `pm::decode_detection_events` on the original detector graph `G`,
  from the full syndrome, on its own instance. Run on **every** shot, not only escalating ones. This
  is the baseline and, inside the model, it is also the core the shot falls back to.
* **`sparse_k`** — the modelled critical path of the decomposed solve on `k` cores:

  ```
  sparse_k = uf + balance + max_c (build_c + solve_c + extract_c) + combine
  ```

  `uf` is the union-find over `H`'s defect-defect edges; `balance` is the cost lookup, the LPT
  assignment and the scatter of defects and edges into the per-core slices; `build_c` is sub-`H_c`
  construction plus `BallMwpm::rebuild` plus the detection-event fill; `solve_c` is
  `process_timeline_until_horizon` alone; `extract_c` is the extraction alone, on a `COMPLETE` core
  only; `combine` is the XOR of the `k` partial observable masks, the sum of the partial weights and
  the OR of the truncation flags. The per-core term is a **max**, never a sum.
* **`system`** — `escalated ? fallback : min(fallback, sparse_k)`. Nothing else in this artifact is
  called `system`. A shot escalates when any core's timeline is `TRUNCATED`, which is the same thing
  as any component truncating, because a core's timeline is the disjoint union of its components'.

`speedup` is `mean(fallback) / mean(system)` over all one million shots of the cell, escalating ones
included.

Syndrome sampling and the ball-graph build are the **input**: they sit outside every timed region.
So does `reset_for_next_shot`, which prepares the instance for the next shot and which this shot's
answer does not depend on — it still runs, in the same place in the same order, immediately after
the extraction region closes.

### The model, and what it does not charge

The `k` loads are solved **serially on one thread** and reported as a max. No thread is ever
spawned. Three things follow, and they all cut the same way:

1. **No migration cost.** Moving a core's slice to that core is not modelled. The scatter happens
   once, serially, inside `balance` on the profiling thread, which is the conservative direction for
   `balance` and the optimistic one for everything downstream.
2. **A warm cache.** A load solved immediately after the previous load sees a warmer cache than a
   real core would. This travels with every number in the file.
3. **A common start.** The fallback core and the `k` solver cores are assumed to start at the same
   instant, when the edge list is ready, and to share nothing.

So `sparse_k` is a floor on what `k` real cores would deliver, not a measurement of them.

### Sweep

| | |
|---|---|
| binary | `load_balancer_profiler`, git `76efe19`, stim v1.16.0, macOS |
| circuit | `surface_code:rotated_memory_x`, `rounds = d`, one observable |
| distances `d` | 17, 21, 25, 27, 29, 31 |
| error rates `p` | 0.0005, 0.001, 0.003 |
| horizon `T` | 1.5, 2.0 (DEM weight units, `T_unit` = median discretised edge of `G`) |
| cores `k` | 5, 9 |
| balancer | LPT, descending component size, cost `s^alpha` with `alpha = 1` |
| shots | 1,000,000 per cell (72 cells, 72,000,000 rows), 1,000 warm-up per cell |
| seed | 20260906 — shot `i` of `(d, p)` is identical across the `T` and `k` sweeps |
| timer | `clock_gettime(CLOCK_THREAD_CPUTIME_ID)`, thread-scoped, 122.944 ns per read |
| correction | none — every row is a raw uncorrected tick delta |
| ran | 2026-09-07 22:54 → 2026-09-08 10:56, `status=ok` |

**No shot is dropped.** The timer is thread-scoped CPU time, so a descheduled shot does not inflate
its own measurement, and there is no contamination probe and no filtered population — every cell
below is the full million. Warm-up shots run the identical workload on the same buffers and read no
clock at all: `run_shot` is templated on the switch, so the reads compile out of the warm-up
instantiation rather than being taken and discarded.

## Results

Speedup ranges from **1.00x** (every `p = 0.003, T = 2` cell) to **5.29x** (`d = 31, p = 0.0005,
T = 1.5, k = 9`). Four things run through the grid.

**It grows with `d`, and only in the sparse regime.** At `p = 0.0005, T = 1.5, k = 9` the speedup
goes 1.65 → 2.35 → 4.34 → 5.00 → 5.08 → 5.29 across `d = 17, 21, 25, 27, 29, 31`. `fallback` scales
with the whole graph; the critical path scales with the largest component `H` breaks the shot into,
and at low `p` that component grows far more slowly than `G` does — mean largest goes 6.2 → 10.8
over the same range of `d` while `fallback` goes 4.4 µs → 57.4 µs.

**It collapses with `p`, and past `p = 0.003` the decomposition is a net loss.** At `p = 0.003,
T = 2` the sparsified graph is not sparse: mean largest component reaches 1,383 defects at `d = 31`
and the critical core holds exactly one component. There `sparse_k` is **1.8x to 2.3x slower than
stock** — 660.9 µs against 374.2 µs at `d = 31, k = 5` — so the per-shot `min` picks `fallback` on
essentially every shot and the speedup column reads a flat 1.00. At `p = 0.003, T = 1.5` the two
sides are a wash (`fallback / sparse_k` between 0.86 and 1.08), and the 1.05–1.12 that shows up in
the speedup column is bought almost entirely by the per-shot `min`: at `d = 25, T = 1.5, k = 5` the
system mean of 136.5 µs is below *both* the fallback mean (142.6) and the sparse mean (156.5),
because on any given shot it takes whichever side was faster.

**`T = 1.5` beats `T = 2` at all 36 `(d, p, k)` points.** The tighter horizon truncates more
aggressively, which is where the win comes from, and it pays for it in escalation: 0.33%–28.8% at
`T = 1.5` against 0.004%–2.56% at `T = 2`. Escalation is cheap in this model — an escalating shot
loses a race it was running anyway, it does not buy a second decode — so the exchange is favourable
everywhere in the sweep. The escalation count is a function of `(d, p, T)` alone and is identical
across `k` by construction; the `k` sweep decodes the same shots the same way and only re-partitions
the work.

**`k = 9` is never worse than `k = 5` on the speedup column, but most of that is not the balancer.**
See the two sections below: the `fallback` control drifted between the two runs, and the balancer's
own contribution — the change in `sparse_k` — is at best 14% and is negative at `T = 2` and at
`p = 0.003`.

### Speedup table

All latencies in microseconds, means over one million shots. `imbalance` is
`k * pred_max / pred_sum` under the LPT prediction, where 1.0 is a perfect split and `k` is
everything on one core. The `x` columns are ratios of `fallback` to `system` at that statistic.

| d | p | T | k | escalated | imbalance | fallback | sparse_k | system | **mean x** | p50 x | p99 x | p99.9 x |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 0.325% | 1.14 | 4.40 | 2.76 | 2.77 | **1.59** | 1.60 | 1.64 | 1.31 |
| 17 | 0.0005 | 1.5 | 9 | 0.325% | 1.54 | 4.47 | 2.70 | 2.71 | **1.65** | 1.67 | 1.64 | 1.30 |
| 17 | 0.0005 | 2 | 5 | 0.004% | 1.36 | 4.47 | 3.54 | 3.42 | **1.31** | 1.32 | 1.28 | 1.36 |
| 17 | 0.0005 | 2 | 9 | 0.004% | 2.25 | 4.58 | 3.58 | 3.45 | **1.33** | 1.36 | 1.28 | 1.36 |
| 17 | 0.001 | 1.5 | 5 | 1.337% | 1.11 | 9.27 | 5.13 | 5.23 | **1.77** | 1.84 | 1.35 | 1.15 |
| 17 | 0.001 | 1.5 | 9 | 1.337% | 1.56 | 9.43 | 4.92 | 5.03 | **1.87** | 1.98 | 1.35 | 1.15 |
| 17 | 0.001 | 2 | 5 | 0.022% | 1.92 | 9.47 | 8.96 | 7.94 | **1.19** | 1.21 | 1.09 | 1.12 |
| 17 | 0.001 | 2 | 9 | 0.022% | 3.41 | 9.67 | 9.19 | 8.12 | **1.19** | 1.20 | 1.09 | 1.12 |
| 17 | 0.003 | 1.5 | 5 | 7.826% | 3.16 | 36.14 | 38.85 | 33.36 | **1.08** | 1.06 | 1.01 | 1.01 |
| 17 | 0.003 | 1.5 | 9 | 7.826% | 5.69 | 36.70 | 39.49 | 33.89 | **1.08** | 1.06 | 1.01 | 1.01 |
| 17 | 0.003 | 2 | 5 | 0.621% | 4.91 | 36.86 | 77.68 | 36.86 | **1.00** | 1.00 | 1.00 | 1.00 |
| 17 | 0.003 | 2 | 9 | 0.621% | 8.83 | 36.94 | 78.09 | 36.94 | **1.00** | 1.00 | 1.00 | 1.00 |
| 21 | 0.0005 | 1.5 | 5 | 0.512% | 1.05 | 8.77 | 4.36 | 4.39 | **2.00** | 2.00 | 1.96 | 1.35 |
| 21 | 0.0005 | 1.5 | 9 | 0.512% | 1.17 | 9.88 | 4.15 | 4.20 | **2.35** | 2.38 | 2.20 | 1.36 |
| 21 | 0.0005 | 2 | 5 | 0.005% | 1.12 | 9.17 | 5.43 | 5.41 | **1.69** | 1.73 | 1.52 | 1.54 |
| 21 | 0.0005 | 2 | 9 | 0.005% | 1.61 | 10.40 | 5.31 | 5.30 | **1.96** | 2.05 | 1.61 | 1.58 |
| 21 | 0.001 | 1.5 | 5 | 2.132% | 1.03 | 18.45 | 8.74 | 9.04 | **2.04** | 2.11 | 1.29 | 1.18 |
| 21 | 0.001 | 1.5 | 9 | 2.132% | 1.20 | 20.27 | 7.91 | 8.27 | **2.45** | 2.61 | 1.30 | 1.18 |
| 21 | 0.001 | 2 | 5 | 0.041% | 1.66 | 19.25 | 15.49 | 14.57 | **1.32** | 1.37 | 1.11 | 1.15 |
| 21 | 0.001 | 2 | 9 | 0.041% | 2.91 | 21.41 | 15.73 | 15.07 | **1.42** | 1.49 | 1.12 | 1.14 |
| 21 | 0.003 | 1.5 | 5 | 12.635% | 3.36 | 72.84 | 84.22 | 69.47 | **1.05** | 1.03 | 1.01 | 1.01 |
| 21 | 0.003 | 1.5 | 9 | 12.635% | 6.05 | 78.14 | 86.44 | 73.72 | **1.06** | 1.04 | 1.01 | 1.01 |
| 21 | 0.003 | 2 | 5 | 1.029% | 4.94 | 74.74 | 169.26 | 74.74 | **1.00** | 1.00 | 1.00 | 1.00 |
| 21 | 0.003 | 2 | 9 | 1.029% | 8.89 | 74.97 | 169.65 | 74.97 | **1.00** | 1.00 | 1.00 | 1.00 |
| 25 | 0.0005 | 1.5 | 5 | 0.768% | 1.02 | 18.14 | 6.66 | 6.78 | **2.67** | 2.71 | 2.26 | 1.30 |
| 25 | 0.0005 | 1.5 | 9 | 0.768% | 1.06 | 27.47 | 6.12 | 6.33 | **4.34** | 4.57 | 3.16 | 1.16 |
| 25 | 0.0005 | 2 | 5 | 0.008% | 1.04 | 19.04 | 8.01 | 8.01 | **2.38** | 2.41 | 2.13 | 1.89 |
| 25 | 0.0005 | 2 | 9 | 0.008% | 1.27 | 27.99 | 7.43 | 7.43 | **3.77** | 3.96 | 2.77 | 2.33 |
| 25 | 0.001 | 1.5 | 5 | 3.114% | 1.01 | 37.45 | 13.97 | 14.85 | **2.52** | 2.67 | 1.22 | 1.16 |
| 25 | 0.001 | 1.5 | 9 | 3.114% | 1.06 | 55.58 | 12.10 | 13.64 | **4.07** | 4.66 | 1.15 | 1.07 |
| 25 | 0.001 | 2 | 5 | 0.057% | 1.51 | 40.22 | 24.31 | 24.05 | **1.67** | 1.83 | 1.16 | 1.17 |
| 25 | 0.001 | 2 | 9 | 0.057% | 2.61 | 56.88 | 24.72 | 24.74 | **2.30** | 2.53 | 1.29 | 1.15 |
| 25 | 0.003 | 1.5 | 5 | 18.420% | 3.59 | 142.63 | 156.52 | 136.49 | **1.05** | 1.02 | 1.01 | 1.01 |
| 25 | 0.003 | 1.5 | 9 | 18.420% | 6.45 | 177.38 | 158.87 | 159.67 | **1.11** | 1.06 | 1.01 | 1.01 |
| 25 | 0.003 | 2 | 5 | 1.573% | 4.96 | 141.44 | 301.28 | 141.44 | **1.00** | 1.00 | 1.00 | 1.00 |
| 25 | 0.003 | 2 | 9 | 1.573% | 8.92 | 142.92 | 302.90 | 142.92 | **1.00** | 1.00 | 1.00 | 1.00 |
| 27 | 0.0005 | 1.5 | 5 | 0.896% | 1.02 | 26.59 | 8.02 | 8.23 | **3.23** | 3.32 | 2.28 | 1.24 |
| 27 | 0.0005 | 1.5 | 9 | 0.896% | 1.04 | 38.80 | 7.42 | 7.76 | **5.00** | 5.32 | 3.02 | 1.14 |
| 27 | 0.0005 | 2 | 5 | 0.008% | 1.03 | 27.35 | 9.60 | 9.60 | **2.85** | 2.88 | 2.56 | 2.20 |
| 27 | 0.0005 | 2 | 9 | 0.008% | 1.18 | 38.84 | 8.86 | 8.86 | **4.38** | 4.58 | 3.24 | 2.76 |
| 27 | 0.001 | 1.5 | 5 | 3.642% | 1.01 | 55.40 | 17.23 | 18.81 | **2.95** | 3.20 | 1.22 | 1.19 |
| 27 | 0.001 | 1.5 | 9 | 3.642% | 1.04 | 77.83 | 14.79 | 17.33 | **4.49** | 5.32 | 1.11 | 1.07 |
| 27 | 0.001 | 2 | 5 | 0.067% | 1.46 | 57.40 | 29.59 | 29.53 | **1.94** | 2.16 | 1.20 | 1.18 |
| 27 | 0.001 | 2 | 9 | 0.067% | 2.50 | 78.54 | 30.48 | 30.53 | **2.57** | 2.82 | 1.41 | 1.19 |
| 27 | 0.003 | 1.5 | 5 | 21.749% | 3.69 | 201.93 | 204.75 | 191.49 | **1.05** | 1.02 | 1.01 | 1.01 |
| 27 | 0.003 | 1.5 | 9 | 21.749% | 6.65 | 243.70 | 210.83 | 218.15 | **1.12** | 1.07 | 1.01 | 1.01 |
| 27 | 0.003 | 2 | 5 | 1.878% | 4.96 | 193.81 | 392.63 | 193.81 | **1.00** | 1.00 | 1.00 | 1.00 |
| 27 | 0.003 | 2 | 9 | 1.878% | 8.93 | 196.34 | 392.12 | 196.34 | **1.00** | 1.00 | 1.00 | 1.00 |
| 29 | 0.0005 | 1.5 | 5 | 1.038% | 1.01 | 40.90 | 9.82 | 10.20 | **4.01** | 4.19 | 1.55 | 1.15 |
| 29 | 0.0005 | 1.5 | 9 | 1.038% | 1.03 | 50.72 | 9.48 | 9.98 | **5.08** | 5.41 | 1.49 | 1.11 |
| 29 | 0.0005 | 2 | 5 | 0.011% | 1.02 | 41.64 | 11.71 | 11.72 | **3.55** | 3.61 | 3.14 | 2.63 |
| 29 | 0.0005 | 2 | 9 | 0.011% | 1.11 | 50.53 | 10.99 | 11.00 | **4.59** | 4.75 | 3.53 | 2.99 |
| 29 | 0.001 | 1.5 | 5 | 4.233% | 1.01 | 83.20 | 21.06 | 23.94 | **3.48** | 3.95 | 1.12 | 1.08 |
| 29 | 0.001 | 1.5 | 9 | 4.233% | 1.02 | 100.73 | 18.22 | 22.01 | **4.58** | 5.57 | 1.09 | 1.05 |
| 29 | 0.001 | 2 | 5 | 0.081% | 1.43 | 86.07 | 36.12 | 36.17 | **2.38** | 2.66 | 1.35 | 1.17 |
| 29 | 0.001 | 2 | 9 | 0.081% | 2.42 | 101.62 | 37.46 | 37.53 | **2.71** | 2.96 | 1.48 | 1.21 |
| 29 | 0.003 | 1.5 | 5 | 25.119% | 3.79 | 280.11 | 265.73 | 263.25 | **1.06** | 1.03 | 1.01 | 1.01 |
| 29 | 0.003 | 1.5 | 9 | 25.119% | 6.82 | 315.54 | 272.36 | 284.61 | **1.11** | 1.06 | 1.01 | 1.01 |
| 29 | 0.003 | 2 | 5 | 2.221% | 4.96 | 264.85 | 502.16 | 264.85 | **1.00** | 1.00 | 1.00 | 1.00 |
| 29 | 0.003 | 2 | 9 | 2.221% | 8.94 | 269.99 | 502.37 | 269.98 | **1.00** | 1.00 | 1.00 | 1.00 |
| 31 | 0.0005 | 1.5 | 5 | 1.207% | 1.01 | 57.38 | 11.93 | 12.55 | **4.57** | 4.86 | 1.32 | 1.12 |
| 31 | 0.0005 | 1.5 | 9 | 1.207% | 1.02 | 63.06 | 11.21 | 11.92 | **5.29** | 5.69 | 1.30 | 1.10 |
| 31 | 0.0005 | 2 | 5 | 0.011% | 1.01 | 57.49 | 14.19 | 14.19 | **4.05** | 4.10 | 3.59 | 3.01 |
| 31 | 0.0005 | 2 | 9 | 0.011% | 1.07 | 63.00 | 12.87 | 12.88 | **4.89** | 5.02 | 3.88 | 3.31 |
| 31 | 0.001 | 1.5 | 5 | 4.893% | 1.00 | 114.51 | 25.74 | 30.40 | **3.77** | 4.46 | 1.09 | 1.06 |
| 31 | 0.001 | 1.5 | 9 | 4.893% | 1.01 | 125.79 | 22.18 | 27.61 | **4.56** | 5.70 | 1.08 | 1.05 |
| 31 | 0.001 | 2 | 5 | 0.094% | 1.39 | 116.94 | 44.35 | 44.43 | **2.63** | 2.94 | 1.49 | 1.20 |
| 31 | 0.001 | 2 | 9 | 0.094% | 2.35 | 125.93 | 45.18 | 45.27 | **2.78** | 3.03 | 1.53 | 1.21 |
| 31 | 0.003 | 1.5 | 5 | 28.783% | 3.88 | 383.39 | 356.08 | 359.46 | **1.07** | 1.03 | 1.15 | 1.00 |
| 31 | 0.003 | 1.5 | 9 | 28.783% | 6.98 | 408.39 | 359.07 | 374.75 | **1.09** | 1.05 | 1.04 | 1.24 |
| 31 | 0.003 | 2 | 5 | 2.561% | 4.97 | 374.19 | 660.88 | 373.83 | **1.00** | 1.00 | 1.00 | 1.03 |
| 31 | 0.003 | 2 | 9 | 2.561% | 8.94 | 363.94 | 647.10 | 363.84 | **1.00** | 1.00 | 1.00 | 1.07 |

### Tail latencies, in microseconds

The mean is not what sizes a syndrome buffer. Blocks arrive on the hardware's clock whether or not
the previous one has been decoded, so the shape of the distribution is what matters, and the tail is
where the two series converge — at `p = 0.003` the p99 ratio is 1.00–1.15 against a mean ratio of
1.00–1.12. That is the expected shape and not a defect: a shot in the tail is a shot with a large
component in `H`, which is exactly the shot that truncates, and on an escalating shot `system` *is*
`fallback` by construction. The decomposed path cannot beat stock on the shots it escalates; it can
only avoid paying for them twice.

The ratio holds up best at low `p`, where escalation is rare enough that even the 99.9th percentile
is still a non-escalating shot: `d = 31, p = 0.0005, T = 2, k = 9` keeps 3.31x at p99.9.

The `max` columns are single worst shots out of a million and are the noisiest thing in the file.
Two of them — `sparse_k` max 4,687.6 µs at `d = 29, p = 0.003, T = 1.5, k = 5` and 4,078.7 µs at
`d = 31` — are an order of magnitude above their own p99.9 and are the profiling thread being
interrupted, not the decoder.

| d | p | T | k | fb p50 | sys p50 | fb p99 | sys p99 | fb p99.9 | sys p99.9 | fb max | sparse_k max |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 4.21 | 2.62 | 8.96 | 5.46 | 12.67 | 9.67 | 45.3 | 25.8 |
| 17 | 0.0005 | 1.5 | 9 | 4.25 | 2.54 | 9.08 | 5.54 | 12.62 | 9.71 | 35.8 | 34.3 |
| 17 | 0.0005 | 2 | 5 | 4.25 | 3.21 | 9.08 | 7.08 | 12.67 | 9.29 | 51.7 | 36.5 |
| 17 | 0.0005 | 2 | 9 | 4.38 | 3.21 | 9.29 | 7.25 | 12.88 | 9.50 | 39.5 | 40.2 |
| 17 | 0.001 | 1.5 | 5 | 8.96 | 4.88 | 16.71 | 12.42 | 22.12 | 19.17 | 54.0 | 46.1 |
| 17 | 0.001 | 1.5 | 9 | 9.08 | 4.58 | 17.04 | 12.62 | 22.42 | 19.42 | 55.0 | 44.2 |
| 17 | 0.001 | 2 | 5 | 9.17 | 7.58 | 17.08 | 15.62 | 22.54 | 20.04 | 58.6 | 56.8 |
| 17 | 0.001 | 2 | 9 | 9.33 | 7.79 | 17.46 | 15.96 | 22.92 | 20.42 | 54.8 | 85.4 |
| 17 | 0.003 | 1.5 | 5 | 35.08 | 33.04 | 59.92 | 59.21 | 76.96 | 76.17 | 257.0 | 150.4 |
| 17 | 0.003 | 1.5 | 9 | 35.67 | 33.58 | 60.71 | 59.96 | 77.96 | 77.25 | 245.2 | 164.5 |
| 17 | 0.003 | 2 | 5 | 35.83 | 35.83 | 60.83 | 60.83 | 77.88 | 77.83 | 242.7 | 241.8 |
| 17 | 0.003 | 2 | 9 | 35.92 | 35.92 | 60.96 | 60.96 | 78.17 | 78.12 | 239.5 | 256.0 |
| 21 | 0.0005 | 1.5 | 5 | 8.50 | 4.25 | 15.58 | 7.96 | 20.96 | 15.50 | 78.8 | 35.7 |
| 21 | 0.0005 | 1.5 | 9 | 9.54 | 4.00 | 17.75 | 8.08 | 23.25 | 17.08 | 192.5 | 35.8 |
| 21 | 0.0005 | 2 | 5 | 8.88 | 5.12 | 16.25 | 10.71 | 21.54 | 13.96 | 70.2 | 46.3 |
| 21 | 0.0005 | 2 | 9 | 10.08 | 4.92 | 18.21 | 11.29 | 23.50 | 14.92 | 76.4 | 50.5 |
| 21 | 0.001 | 1.5 | 5 | 18.04 | 8.54 | 29.42 | 22.88 | 37.33 | 31.62 | 90.3 | 43.2 |
| 21 | 0.001 | 1.5 | 9 | 19.79 | 7.58 | 32.25 | 24.79 | 40.25 | 34.21 | 125.9 | 48.5 |
| 21 | 0.001 | 2 | 5 | 18.83 | 13.71 | 30.50 | 27.58 | 38.12 | 33.29 | 95.2 | 87.9 |
| 21 | 0.001 | 2 | 9 | 20.96 | 14.08 | 33.58 | 29.92 | 41.04 | 36.08 | 219.9 | 83.6 |
| 21 | 0.003 | 1.5 | 5 | 71.58 | 69.58 | 106.67 | 105.92 | 128.17 | 126.62 | 615.4 | 255.2 |
| 21 | 0.003 | 1.5 | 9 | 76.83 | 74.17 | 113.58 | 112.62 | 135.96 | 134.25 | 625.4 | 330.8 |
| 21 | 0.003 | 2 | 5 | 73.46 | 73.46 | 108.75 | 108.75 | 130.29 | 130.29 | 662.9 | 632.1 |
| 21 | 0.003 | 2 | 9 | 73.71 | 73.71 | 108.96 | 108.96 | 130.50 | 130.50 | 623.5 | 456.5 |
| 25 | 0.0005 | 1.5 | 5 | 17.75 | 6.54 | 28.83 | 12.75 | 35.67 | 27.54 | 86.9 | 38.7 |
| 25 | 0.0005 | 1.5 | 9 | 27.21 | 5.96 | 39.62 | 12.54 | 45.71 | 39.33 | 114.0 | 67.6 |
| 25 | 0.0005 | 2 | 5 | 18.67 | 7.75 | 29.79 | 14.00 | 36.29 | 19.25 | 84.4 | 69.1 |
| 25 | 0.0005 | 2 | 9 | 27.71 | 7.00 | 40.21 | 14.50 | 46.25 | 19.87 | 88.5 | 45.5 |
| 25 | 0.001 | 1.5 | 5 | 36.88 | 13.83 | 54.54 | 44.54 | 65.67 | 56.38 | 174.7 | 74.4 |
| 25 | 0.001 | 1.5 | 9 | 55.17 | 11.83 | 74.46 | 64.79 | 83.54 | 77.83 | 210.7 | 64.5 |
| 25 | 0.001 | 2 | 5 | 39.67 | 21.62 | 57.58 | 49.71 | 68.58 | 58.54 | 169.7 | 289.0 |
| 25 | 0.001 | 2 | 9 | 56.46 | 22.29 | 76.17 | 59.12 | 85.54 | 74.17 | 370.1 | 137.0 |
| 25 | 0.003 | 1.5 | 5 | 140.92 | 137.83 | 192.25 | 190.71 | 223.54 | 220.42 | 611.5 | 679.5 |
| 25 | 0.003 | 1.5 | 9 | 175.96 | 166.12 | 228.67 | 226.67 | 257.42 | 254.67 | 704.3 | 566.3 |
| 25 | 0.003 | 2 | 5 | 139.75 | 139.75 | 190.04 | 190.04 | 220.79 | 220.79 | 621.5 | 913.3 |
| 25 | 0.003 | 2 | 9 | 141.25 | 141.25 | 191.29 | 191.29 | 221.54 | 221.54 | 647.6 | 763.3 |
| 27 | 0.0005 | 1.5 | 5 | 26.17 | 7.88 | 39.62 | 17.42 | 47.08 | 37.83 | 102.7 | 72.4 |
| 27 | 0.0005 | 1.5 | 9 | 38.54 | 7.25 | 53.29 | 17.67 | 60.12 | 52.88 | 130.7 | 84.7 |
| 27 | 0.0005 | 2 | 5 | 26.96 | 9.37 | 40.17 | 15.67 | 47.58 | 21.62 | 126.2 | 109.4 |
| 27 | 0.0005 | 2 | 9 | 38.54 | 8.42 | 53.46 | 16.50 | 61.62 | 22.29 | 224.9 | 56.6 |
| 27 | 0.001 | 1.5 | 5 | 54.50 | 17.04 | 78.92 | 64.71 | 95.42 | 80.17 | 448.9 | 187.2 |
| 27 | 0.001 | 1.5 | 9 | 77.42 | 14.54 | 100.04 | 89.75 | 110.58 | 103.79 | 249.5 | 91.9 |
| 27 | 0.001 | 2 | 5 | 56.79 | 26.33 | 78.42 | 65.08 | 90.58 | 76.88 | 214.8 | 320.1 |
| 27 | 0.001 | 2 | 9 | 78.12 | 27.67 | 100.88 | 71.38 | 110.83 | 93.21 | 526.8 | 428.4 |
| 27 | 0.003 | 1.5 | 5 | 200.12 | 195.67 | 260.00 | 258.08 | 294.58 | 290.50 | 703.5 | 785.8 |
| 27 | 0.003 | 1.5 | 9 | 242.25 | 227.25 | 303.04 | 300.50 | 334.50 | 331.58 | 725.8 | 687.1 |
| 27 | 0.003 | 2 | 5 | 191.96 | 191.96 | 251.25 | 251.25 | 286.33 | 286.33 | 707.2 | 1000.7 |
| 27 | 0.003 | 2 | 9 | 194.50 | 194.50 | 253.12 | 253.12 | 286.96 | 286.96 | 664.2 | 948.7 |
| 29 | 0.0005 | 1.5 | 5 | 40.54 | 9.67 | 56.12 | 36.29 | 63.62 | 55.29 | 146.7 | 129.1 |
| 29 | 0.0005 | 1.5 | 9 | 50.46 | 9.33 | 67.38 | 45.29 | 74.54 | 67.00 | 169.3 | 71.4 |
| 29 | 0.0005 | 2 | 5 | 41.33 | 11.46 | 56.83 | 18.08 | 64.25 | 24.46 | 470.6 | 102.6 |
| 29 | 0.0005 | 2 | 9 | 50.25 | 10.58 | 67.12 | 19.04 | 74.42 | 24.92 | 331.0 | 61.2 |
| 29 | 0.001 | 1.5 | 5 | 82.67 | 20.92 | 107.38 | 95.58 | 119.54 | 110.50 | 303.6 | 190.1 |
| 29 | 0.001 | 1.5 | 9 | 100.33 | 18.00 | 125.71 | 115.08 | 136.67 | 130.25 | 287.8 | 117.0 |
| 29 | 0.001 | 2 | 5 | 85.58 | 32.21 | 110.25 | 81.37 | 122.04 | 104.04 | 416.4 | 514.0 |
| 29 | 0.001 | 2 | 9 | 101.21 | 34.17 | 127.00 | 85.54 | 138.04 | 113.75 | 298.1 | 282.9 |
| 29 | 0.003 | 1.5 | 5 | 278.42 | 271.29 | 346.29 | 344.29 | 383.46 | 379.71 | 746.6 | 4687.6 |
| 29 | 0.003 | 1.5 | 9 | 314.00 | 295.83 | 381.83 | 379.04 | 416.25 | 412.67 | 825.9 | 789.6 |
| 29 | 0.003 | 2 | 5 | 262.96 | 262.96 | 330.71 | 330.71 | 369.58 | 369.58 | 835.2 | 1405.0 |
| 29 | 0.003 | 2 | 9 | 268.17 | 268.17 | 335.38 | 335.38 | 372.42 | 372.42 | 950.1 | 1462.3 |
| 31 | 0.0005 | 1.5 | 5 | 57.08 | 11.75 | 74.79 | 56.67 | 82.71 | 74.12 | 193.3 | 144.1 |
| 31 | 0.0005 | 1.5 | 9 | 62.79 | 11.04 | 81.50 | 62.71 | 89.46 | 81.12 | 376.4 | 54.3 |
| 31 | 0.0005 | 2 | 5 | 57.17 | 13.96 | 74.96 | 20.87 | 83.00 | 27.58 | 322.0 | 168.1 |
| 31 | 0.0005 | 2 | 9 | 62.71 | 12.50 | 81.38 | 21.00 | 89.42 | 27.04 | 286.3 | 132.2 |
| 31 | 0.001 | 1.5 | 5 | 114.08 | 25.58 | 141.46 | 129.83 | 153.83 | 145.62 | 553.4 | 100.5 |
| 31 | 0.001 | 1.5 | 9 | 125.33 | 22.00 | 153.92 | 142.58 | 166.38 | 158.96 | 575.7 | 142.0 |
| 31 | 0.001 | 2 | 5 | 116.50 | 39.58 | 144.54 | 97.17 | 156.62 | 130.58 | 572.0 | 464.5 |
| 31 | 0.001 | 2 | 9 | 125.50 | 41.38 | 153.96 | 100.83 | 166.42 | 137.00 | 596.3 | 203.5 |
| 31 | 0.003 | 1.5 | 5 | 370.58 | 360.00 | 626.92 | 542.92 | 1745.96 | 1738.33 | 3519.8 | 4078.7 |
| 31 | 0.003 | 1.5 | 9 | 404.50 | 386.21 | 514.67 | 495.71 | 891.04 | 716.21 | 1768.9 | 1558.5 |
| 31 | 0.003 | 2 | 5 | 363.75 | 363.75 | 595.12 | 593.38 | 1165.75 | 1129.50 | 2180.8 | 2723.4 |
| 31 | 0.003 | 2 | 9 | 358.83 | 358.83 | 487.71 | 487.71 | 805.29 | 749.54 | 1949.9 | 2814.8 |

### Where the critical path goes

Module means in raw nanoseconds, from the 128-bit accumulators in `agg.json`. The six sum to the
`sparse_k` mean by construction. `serial share` is `(uf + balance + combine) / sparse_k` — the part
that does not parallelise no matter how many cores are added.

Two readings dominate:

**`build_crit` is the largest single module wherever the components are big.** At `T = 2` it exceeds
`solve_crit` in every cell, by up to 1.86x, and at `p = 0.003, T = 2` it is 50–58% of the whole
critical path. Sub-`H_c` construction and `BallMwpm::rebuild` cost more than the blossom solve they
set up. If anything on this path is worth optimising, it is that.

**The serial prologue is what caps the `k` sweep.** `balance` grows with `k` — LPT scans the `k`
accumulated costs per component — while the parallel part shrinks, so the serial share climbs from
7.5% (`d = 31, p = 0.003, T = 2, k = 5`) to 51.7% (`d = 31, p = 0.0005, T = 1.5, k = 9`), where
`balance` alone is 32.8% of the path. That is the whole Amdahl story of the `k` column: at low `p`,
nine cores spend more than half the critical path deciding what the nine cores should do.

`combine` is 127–172 ns everywhere, independent of `d`, `p`, `T` and `k`, and is under 5% of the
path in every cell.

| d | p | T | k | uf | balance | build_crit | solve_crit | extract_crit | combine | sum = sparse_k | serial share |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 372 | 511 | 705 | 816 | 233 | 127 | 2,764 | 36.5% |
| 17 | 0.0005 | 1.5 | 9 | 375 | 665 | 633 | 682 | 213 | 135 | 2,704 | 43.5% |
| 17 | 0.0005 | 2 | 5 | 424 | 538 | 1,147 | 1,049 | 255 | 128 | 3,542 | 30.8% |
| 17 | 0.0005 | 2 | 9 | 427 | 659 | 1,100 | 1,000 | 258 | 136 | 3,578 | 34.1% |
| 17 | 0.001 | 1.5 | 5 | 709 | 895 | 1,369 | 1,679 | 344 | 130 | 5,126 | 33.8% |
| 17 | 0.001 | 1.5 | 9 | 713 | 1,102 | 1,253 | 1,418 | 300 | 138 | 4,924 | 39.7% |
| 17 | 0.001 | 2 | 5 | 829 | 1,067 | 3,415 | 3,011 | 504 | 130 | 8,957 | 22.6% |
| 17 | 0.001 | 2 | 9 | 836 | 1,190 | 3,521 | 3,000 | 505 | 138 | 9,190 | 23.5% |
| 17 | 0.003 | 1.5 | 5 | 2,431 | 3,209 | 15,271 | 15,873 | 1,929 | 132 | 38,845 | 14.9% |
| 17 | 0.003 | 1.5 | 9 | 2,456 | 3,334 | 15,644 | 15,985 | 1,931 | 140 | 39,491 | 15.0% |
| 17 | 0.003 | 2 | 5 | 3,009 | 4,252 | 39,165 | 28,035 | 3,086 | 128 | 77,676 | 9.5% |
| 17 | 0.003 | 2 | 9 | 3,013 | 4,293 | 39,455 | 28,122 | 3,080 | 131 | 78,094 | 9.5% |
| 21 | 0.0005 | 1.5 | 5 | 625 | 900 | 1,040 | 1,327 | 319 | 144 | 4,355 | 38.3% |
| 21 | 0.0005 | 1.5 | 9 | 638 | 1,212 | 889 | 1,002 | 256 | 152 | 4,149 | 48.2% |
| 21 | 0.0005 | 2 | 5 | 722 | 945 | 1,694 | 1,588 | 337 | 146 | 5,432 | 33.4% |
| 21 | 0.0005 | 2 | 9 | 731 | 1,151 | 1,588 | 1,383 | 302 | 152 | 5,308 | 38.3% |
| 21 | 0.001 | 1.5 | 5 | 1,285 | 1,648 | 2,387 | 2,775 | 510 | 137 | 8,742 | 35.1% |
| 21 | 0.001 | 1.5 | 9 | 1,305 | 2,020 | 1,984 | 2,053 | 390 | 157 | 7,908 | 44.0% |
| 21 | 0.001 | 2 | 5 | 1,525 | 1,938 | 6,317 | 4,814 | 756 | 137 | 15,488 | 23.2% |
| 21 | 0.001 | 2 | 9 | 1,536 | 2,083 | 6,442 | 4,776 | 740 | 156 | 15,732 | 24.0% |
| 21 | 0.003 | 1.5 | 5 | 4,615 | 6,261 | 36,737 | 32,747 | 3,730 | 135 | 84,225 | 13.1% |
| 21 | 0.003 | 1.5 | 9 | 4,652 | 6,426 | 38,519 | 32,943 | 3,738 | 161 | 86,439 | 13.0% |
| 21 | 0.003 | 2 | 5 | 5,881 | 7,874 | 94,163 | 55,301 | 5,906 | 132 | 169,258 | 8.2% |
| 21 | 0.003 | 2 | 9 | 5,887 | 7,907 | 94,514 | 55,296 | 5,910 | 135 | 169,649 | 8.2% |
| 25 | 0.0005 | 1.5 | 5 | 992 | 1,473 | 1,563 | 2,048 | 447 | 136 | 6,659 | 39.1% |
| 25 | 0.0005 | 1.5 | 9 | 1,009 | 1,960 | 1,189 | 1,491 | 326 | 145 | 6,120 | 50.9% |
| 25 | 0.0005 | 2 | 5 | 1,137 | 1,504 | 2,453 | 2,334 | 449 | 135 | 8,011 | 34.6% |
| 25 | 0.0005 | 2 | 9 | 1,161 | 1,814 | 2,094 | 1,857 | 360 | 144 | 7,430 | 42.0% |
| 25 | 0.001 | 1.5 | 5 | 2,140 | 2,730 | 3,872 | 4,331 | 757 | 135 | 13,966 | 35.8% |
| 25 | 0.001 | 1.5 | 9 | 2,227 | 3,364 | 2,793 | 3,052 | 521 | 147 | 12,104 | 47.4% |
| 25 | 0.001 | 2 | 5 | 2,580 | 3,222 | 9,968 | 7,305 | 1,097 | 136 | 24,308 | 24.4% |
| 25 | 0.001 | 2 | 9 | 2,689 | 3,436 | 10,182 | 7,203 | 1,065 | 147 | 24,721 | 25.4% |
| 25 | 0.003 | 1.5 | 5 | 7,818 | 10,710 | 71,464 | 60,029 | 6,368 | 135 | 156,524 | 11.9% |
| 25 | 0.003 | 1.5 | 9 | 7,803 | 10,908 | 73,392 | 60,285 | 6,336 | 145 | 158,869 | 11.9% |
| 25 | 0.003 | 2 | 5 | 10,052 | 12,687 | 172,854 | 95,563 | 9,989 | 132 | 301,276 | 7.6% |
| 25 | 0.003 | 2 | 9 | 10,104 | 12,795 | 173,810 | 96,009 | 10,046 | 136 | 302,900 | 7.6% |
| 27 | 0.0005 | 1.5 | 5 | 1,205 | 1,817 | 1,873 | 2,473 | 518 | 134 | 8,020 | 39.4% |
| 27 | 0.0005 | 1.5 | 9 | 1,244 | 2,397 | 1,406 | 1,854 | 373 | 146 | 7,421 | 51.0% |
| 27 | 0.0005 | 2 | 5 | 1,383 | 1,860 | 2,922 | 2,781 | 517 | 133 | 9,596 | 35.2% |
| 27 | 0.0005 | 2 | 9 | 1,442 | 2,229 | 2,429 | 2,215 | 397 | 145 | 8,858 | 43.1% |
| 27 | 0.001 | 1.5 | 5 | 2,654 | 3,403 | 4,832 | 5,294 | 913 | 133 | 17,229 | 35.9% |
| 27 | 0.001 | 1.5 | 9 | 2,734 | 4,168 | 3,366 | 3,772 | 599 | 149 | 14,788 | 47.7% |
| 27 | 0.001 | 2 | 5 | 3,203 | 3,987 | 12,211 | 8,754 | 1,299 | 134 | 29,587 | 24.8% |
| 27 | 0.001 | 2 | 9 | 3,266 | 4,233 | 12,596 | 8,997 | 1,242 | 150 | 30,484 | 25.1% |
| 27 | 0.003 | 1.5 | 5 | 9,915 | 13,487 | 94,713 | 78,439 | 8,060 | 136 | 204,749 | 11.5% |
| 27 | 0.003 | 1.5 | 9 | 10,173 | 13,984 | 98,330 | 80,055 | 8,138 | 146 | 210,827 | 11.5% |
| 27 | 0.003 | 2 | 5 | 12,844 | 16,711 | 226,151 | 123,806 | 12,983 | 133 | 392,628 | 7.6% |
| 27 | 0.003 | 2 | 9 | 12,794 | 16,718 | 226,122 | 123,445 | 12,901 | 136 | 392,117 | 7.6% |
| 29 | 0.0005 | 1.5 | 5 | 1,500 | 2,243 | 2,297 | 3,033 | 613 | 136 | 9,821 | 39.5% |
| 29 | 0.0005 | 1.5 | 9 | 1,628 | 3,067 | 1,784 | 2,409 | 428 | 162 | 9,478 | 51.2% |
| 29 | 0.0005 | 2 | 5 | 1,723 | 2,299 | 3,551 | 3,396 | 608 | 135 | 11,711 | 35.5% |
| 29 | 0.0005 | 2 | 9 | 1,857 | 2,863 | 2,966 | 2,702 | 444 | 160 | 10,993 | 44.4% |
| 29 | 0.001 | 1.5 | 5 | 3,228 | 4,202 | 5,892 | 6,497 | 1,101 | 135 | 21,055 | 35.9% |
| 29 | 0.001 | 1.5 | 9 | 3,376 | 5,245 | 4,090 | 4,657 | 696 | 154 | 18,219 | 48.2% |
| 29 | 0.001 | 2 | 5 | 3,917 | 4,967 | 14,927 | 10,632 | 1,545 | 136 | 36,123 | 25.0% |
| 29 | 0.001 | 2 | 9 | 4,079 | 5,370 | 15,505 | 10,875 | 1,476 | 156 | 37,460 | 25.6% |
| 29 | 0.003 | 1.5 | 5 | 12,378 | 17,093 | 124,257 | 101,779 | 10,091 | 134 | 265,733 | 11.1% |
| 29 | 0.003 | 1.5 | 9 | 12,683 | 17,635 | 128,170 | 103,534 | 10,174 | 161 | 272,357 | 11.2% |
| 29 | 0.003 | 2 | 5 | 15,953 | 21,858 | 290,436 | 157,389 | 16,393 | 134 | 502,163 | 7.6% |
| 29 | 0.003 | 2 | 9 | 15,917 | 21,902 | 290,955 | 157,113 | 16,344 | 138 | 502,368 | 7.6% |
| 31 | 0.0005 | 1.5 | 5 | 1,841 | 2,753 | 2,797 | 3,683 | 720 | 135 | 11,928 | 39.6% |
| 31 | 0.0005 | 1.5 | 9 | 1,952 | 3,677 | 2,104 | 2,819 | 496 | 165 | 11,213 | 51.7% |
| 31 | 0.0005 | 2 | 5 | 2,115 | 2,807 | 4,290 | 4,129 | 714 | 133 | 14,188 | 35.6% |
| 31 | 0.0005 | 2 | 9 | 2,235 | 3,416 | 3,413 | 3,138 | 505 | 164 | 12,871 | 45.2% |
| 31 | 0.001 | 1.5 | 5 | 3,970 | 5,135 | 7,287 | 7,890 | 1,326 | 136 | 25,743 | 35.9% |
| 31 | 0.001 | 1.5 | 9 | 4,134 | 6,412 | 5,050 | 5,592 | 821 | 168 | 22,178 | 48.3% |
| 31 | 0.001 | 2 | 5 | 4,849 | 6,109 | 18,412 | 13,005 | 1,837 | 138 | 44,350 | 25.0% |
| 31 | 0.001 | 2 | 9 | 4,951 | 6,489 | 18,999 | 12,838 | 1,732 | 166 | 45,175 | 25.7% |
| 31 | 0.003 | 1.5 | 5 | 16,169 | 22,095 | 167,119 | 137,458 | 13,090 | 146 | 356,078 | 10.8% |
| 31 | 0.003 | 1.5 | 9 | 16,227 | 22,422 | 170,220 | 137,007 | 13,021 | 172 | 359,070 | 10.8% |
| 31 | 0.003 | 2 | 5 | 20,285 | 28,812 | 383,145 | 207,084 | 21,415 | 143 | 660,884 | 7.5% |
| 31 | 0.003 | 2 | 9 | 19,964 | 28,597 | 375,320 | 202,326 | 20,752 | 144 | 647,102 | 7.5% |

### Structure of H, and why the balancer runs out of room

`largest` is the mean largest component of `H`; it and the escalation rate depend on `(d, p, T)`
only. `crit comps` is the mean number of components on the critical core.

LPT cannot split a component, so no assignment can do better than
`k * cost(largest) / total cost`, whatever `k` is. That floor is the whole pattern here: where components are small and numerous
the split is near-perfect (1.00–1.06 at `p ≤ 0.001, T = 1.5, d ≥ 25`), and where one component
dominates the shot the imbalance saturates at `k` — 8.94 out of 9 at `p = 0.003, T = 2`, meaning
99% of the predicted cost lands on one core, which holds exactly 1.00 components. Adding cores
there does nothing, which is exactly what the `sparse_k` column shows.

| d | p | T | escalated | mean largest | imbal k=5 | imbal k=9 | crit comps k=5 | crit comps k=9 |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 0.325% | 6.2 | 1.14 | 1.54 | 2.39 | 1.22 |
| 17 | 0.0005 | 2 | 0.004% | 9.7 | 1.36 | 2.25 | 1.49 | 1.02 |
| 17 | 0.001 | 1.5 | 1.337% | 13.0 | 1.11 | 1.56 | 3.21 | 1.31 |
| 17 | 0.001 | 2 | 0.022% | 29.9 | 1.92 | 3.41 | 1.20 | 1.00 |
| 17 | 0.003 | 1.5 | 7.826% | 142.6 | 3.16 | 5.69 | 1.03 | 1.00 |
| 17 | 0.003 | 2 | 0.621% | 217.5 | 4.91 | 8.83 | 1.00 | 1.00 |
| 21 | 0.0005 | 1.5 | 0.512% | 7.6 | 1.05 | 1.17 | 4.97 | 2.13 |
| 21 | 0.0005 | 2 | 0.005% | 13.0 | 1.12 | 1.61 | 2.88 | 1.22 |
| 21 | 0.001 | 1.5 | 2.132% | 17.2 | 1.03 | 1.20 | 7.03 | 2.60 |
| 21 | 0.001 | 2 | 0.041% | 48.5 | 1.66 | 2.91 | 1.65 | 1.03 |
| 21 | 0.003 | 1.5 | 12.635% | 287.9 | 3.36 | 6.05 | 1.01 | 1.00 |
| 21 | 0.003 | 2 | 1.029% | 418.7 | 4.94 | 8.89 | 1.00 | 1.00 |
| 25 | 0.0005 | 1.5 | 0.768% | 9.0 | 1.02 | 1.06 | 9.00 | 4.29 |
| 25 | 0.0005 | 2 | 0.008% | 16.3 | 1.04 | 1.27 | 5.64 | 2.05 |
| 25 | 0.001 | 1.5 | 3.114% | 21.4 | 1.01 | 1.06 | 12.74 | 5.54 |
| 25 | 0.001 | 2 | 0.057% | 73.5 | 1.51 | 2.61 | 2.46 | 1.08 |
| 25 | 0.003 | 1.5 | 18.420% | 521.3 | 3.59 | 6.45 | 1.00 | 1.00 |
| 25 | 0.003 | 2 | 1.573% | 715.9 | 4.96 | 8.92 | 1.00 | 1.00 |
| 27 | 0.0005 | 1.5 | 0.896% | 9.6 | 1.02 | 1.04 | 11.51 | 5.77 |
| 27 | 0.0005 | 2 | 0.008% | 18.0 | 1.03 | 1.18 | 7.47 | 2.79 |
| 27 | 0.001 | 1.5 | 3.642% | 23.4 | 1.01 | 1.04 | 16.27 | 7.63 |
| 27 | 0.001 | 2 | 0.067% | 88.8 | 1.46 | 2.50 | 2.97 | 1.13 |
| 27 | 0.003 | 1.5 | 21.749% | 677.9 | 3.69 | 6.65 | 1.00 | 1.00 |
| 27 | 0.003 | 2 | 1.878% | 906.3 | 4.96 | 8.93 | 1.00 | 1.00 |
| 29 | 0.0005 | 1.5 | 1.038% | 10.2 | 1.01 | 1.03 | 14.41 | 7.55 |
| 29 | 0.0005 | 2 | 0.011% | 19.7 | 1.02 | 1.11 | 9.60 | 3.87 |
| 29 | 0.001 | 1.5 | 4.233% | 25.5 | 1.01 | 1.02 | 20.24 | 10.04 |
| 29 | 0.001 | 2 | 0.081% | 106.6 | 1.43 | 2.42 | 3.54 | 1.18 |
| 29 | 0.003 | 1.5 | 25.119% | 864.3 | 3.79 | 6.82 | 1.00 | 1.00 |
| 29 | 0.003 | 2 | 2.221% | 1127.8 | 4.96 | 8.94 | 1.00 | 1.00 |
| 31 | 0.0005 | 1.5 | 1.207% | 10.8 | 1.01 | 1.02 | 17.73 | 9.45 |
| 31 | 0.0005 | 2 | 0.011% | 21.4 | 1.01 | 1.07 | 12.03 | 5.15 |
| 31 | 0.001 | 1.5 | 4.893% | 27.5 | 1.00 | 1.01 | 24.77 | 12.79 |
| 31 | 0.001 | 2 | 0.094% | 126.7 | 1.39 | 2.35 | 4.18 | 1.24 |
| 31 | 0.003 | 1.5 | 28.783% | 1083.2 | 3.88 | 6.98 | 1.00 | 1.00 |
| 31 | 0.003 | 2 | 2.561% | 1382.6 | 4.97 | 8.94 | 1.00 | 1.00 |

## Caveats

### The `fallback` control moved with `k`

`fallback` is the internal control of the `k` sweep: it is the same shots, decoded the same way, on
the same instance, in every `k` cell. It should not move with `k`. It does — monotonically upward in
35 of 36 pairs, by up to **1.52x** (`d = 25, p = 0.0005, T = 1.5`: 18.14 µs → 27.47 µs). The `k = 9`
cell of a pair always runs immediately after the `k = 5` cell in the same process, with more solver
instances resident, so this is an ordering or memory-pressure effect on the profiling machine rather
than anything the balancer did.

The consequence is direct: **the `k = 9` speedup column overstates the balancer.** The balancer's
own effect is the last column below — the change in `sparse_k`, which does not depend on the
fallback. It is a gain of at most 14% (`d = 27, p = 0.001, T = 1.5`: 0.858x) and is a *loss* at
every `T = 2` cell and every `p = 0.003` cell. Read the `k` sweep from that column, not from
`speedup`.

| d | p | T | fallback k=5 | fallback k=9 | drift | sparse_k k=5 | sparse_k k=9 | balancer |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 4.40 | 4.47 | 1.016x | 2.76 | 2.70 | 0.978x |
| 17 | 0.0005 | 2 | 4.47 | 4.58 | 1.024x | 3.54 | 3.58 | 1.010x |
| 17 | 0.001 | 1.5 | 9.27 | 9.43 | 1.017x | 5.13 | 4.92 | 0.961x |
| 17 | 0.001 | 2 | 9.47 | 9.67 | 1.020x | 8.96 | 9.19 | 1.026x |
| 17 | 0.003 | 1.5 | 36.14 | 36.70 | 1.016x | 38.85 | 39.49 | 1.017x |
| 17 | 0.003 | 2 | 36.86 | 36.94 | 1.002x | 77.68 | 78.09 | 1.005x |
| 21 | 0.0005 | 1.5 | 8.77 | 9.88 | 1.127x | 4.36 | 4.15 | 0.953x |
| 21 | 0.0005 | 2 | 9.17 | 10.40 | 1.134x | 5.43 | 5.31 | 0.977x |
| 21 | 0.001 | 1.5 | 18.45 | 20.27 | 1.099x | 8.74 | 7.91 | 0.905x |
| 21 | 0.001 | 2 | 19.25 | 21.41 | 1.112x | 15.49 | 15.73 | 1.016x |
| 21 | 0.003 | 1.5 | 72.84 | 78.14 | 1.073x | 84.22 | 86.44 | 1.026x |
| 21 | 0.003 | 2 | 74.74 | 74.97 | 1.003x | 169.26 | 169.65 | 1.002x |
| 25 | 0.0005 | 1.5 | 18.14 | 27.47 | 1.515x | 6.66 | 6.12 | 0.919x |
| 25 | 0.0005 | 2 | 19.04 | 27.99 | 1.470x | 8.01 | 7.43 | 0.927x |
| 25 | 0.001 | 1.5 | 37.45 | 55.58 | 1.484x | 13.97 | 12.10 | 0.867x |
| 25 | 0.001 | 2 | 40.22 | 56.88 | 1.414x | 24.31 | 24.72 | 1.017x |
| 25 | 0.003 | 1.5 | 142.63 | 177.38 | 1.244x | 156.52 | 158.87 | 1.015x |
| 25 | 0.003 | 2 | 141.44 | 142.92 | 1.010x | 301.28 | 302.90 | 1.005x |
| 27 | 0.0005 | 1.5 | 26.59 | 38.80 | 1.459x | 8.02 | 7.42 | 0.925x |
| 27 | 0.0005 | 2 | 27.35 | 38.84 | 1.420x | 9.60 | 8.86 | 0.923x |
| 27 | 0.001 | 1.5 | 55.40 | 77.83 | 1.405x | 17.23 | 14.79 | 0.858x |
| 27 | 0.001 | 2 | 57.40 | 78.54 | 1.368x | 29.59 | 30.48 | 1.030x |
| 27 | 0.003 | 1.5 | 201.93 | 243.70 | 1.207x | 204.75 | 210.83 | 1.030x |
| 27 | 0.003 | 2 | 193.81 | 196.34 | 1.013x | 392.63 | 392.12 | 0.999x |
| 29 | 0.0005 | 1.5 | 40.90 | 50.72 | 1.240x | 9.82 | 9.48 | 0.965x |
| 29 | 0.0005 | 2 | 41.64 | 50.53 | 1.213x | 11.71 | 10.99 | 0.939x |
| 29 | 0.001 | 1.5 | 83.20 | 100.73 | 1.211x | 21.06 | 18.22 | 0.865x |
| 29 | 0.001 | 2 | 86.07 | 101.62 | 1.181x | 36.12 | 37.46 | 1.037x |
| 29 | 0.003 | 1.5 | 280.11 | 315.54 | 1.126x | 265.73 | 272.36 | 1.025x |
| 29 | 0.003 | 2 | 264.85 | 269.99 | 1.019x | 502.16 | 502.37 | 1.000x |
| 31 | 0.0005 | 1.5 | 57.38 | 63.06 | 1.099x | 11.93 | 11.21 | 0.940x |
| 31 | 0.0005 | 2 | 57.49 | 63.00 | 1.096x | 14.19 | 12.87 | 0.907x |
| 31 | 0.001 | 1.5 | 114.51 | 125.79 | 1.099x | 25.74 | 22.18 | 0.861x |
| 31 | 0.001 | 2 | 116.94 | 125.93 | 1.077x | 44.35 | 45.18 | 1.019x |
| 31 | 0.003 | 1.5 | 383.39 | 408.39 | 1.065x | 356.08 | 359.07 | 1.008x |
| 31 | 0.003 | 2 | 374.19 | 363.94 | 0.973x | 660.88 | 647.10 | 0.979x |

### The rows are uncorrected, and the two series are not read-symmetric

Every tick column is a raw delta of a clock that costs **122.944 ns per read**, and a bracketed
region costs about one read. `sparse_k` sums six intervals and therefore carries about six of those;
`fallback` is one interval and carries one. Subtracting a flat overhead from both would leave five
of them inside `sparse_k` and bias the ratio against the decomposed path.

At large `d` and high `p` this is noise. At the small, fast cells it is not: at `d = 17,
p = 0.0005, T = 1.5, k = 5` the six reads are ~738 ns against a `sparse_k` mean of 2,764 ns — **27%
of the measurement**. Correcting per column moves `mean(fallback) / mean(sparse_k)` at that cell
from 1.59 to 2.11; at `d = 31, p = 0.0005, T = 1.5, k = 9` (6.6% overhead) from 5.62 to 6.01. The
plotter does this per column behind `--subtract-overhead`, with the interval count per column and
five rather than six on a shot whose critical core truncated and so had no extraction. It is off
here, and the tables above are the raw numbers.

### Other things travelling with every row

* Single machine, single thread, macOS, `-O3`. The earlier campaign in this project ran on a Linux
  `rdpmc` backend where the timer correction was tens of nanoseconds rather than 123.
* `alpha = 1`. The design specifies a default of 1.1 fitted to per-component solve time; this run
  used a linear cost model on instruction. The balancer's *prediction* is what changes with `alpha`,
  not the measured ticks it is predicting.
* `--verify` and `--check-determinism` were both off for this campaign (`verify=0`,
  `check_determinism=0` in `run.log`).
* `no_growth_after_warmup` is 0 in 29 of the 72 cells — those cells grew a capacity after warm-up,
  so a handful of their early measured shots include an allocation. `run.log` records the flag per
  cell rather than asserting the condition away.

## Reproducing

```
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --target load_balancer_profiler -j4
build/load_balancer_profiler --d 17,21,25,27,29,31 --p 0.0005,0.001,0.003 --T 1.5,2 --k 5,9 \
    --alpha 1 --shots 1000000 --warmup 1000 --seed 20260906 --out results
```

That writes `shots.csv` (one row per shot — 8.2 GB for this sweep), `agg.json` (per-cell means,
percentiles, 128-bit module accumulators and the `largest` / `crit_ncomps` histograms) and `run.log`
(git hash, build flags, stim version, timer backend and overhead, the exact generator call per
corpus, and every region definition quoted in this file). Add `--verify` on a shorter run to check
the decomposition against the monolithic per-component path shot by shot, and `--check-determinism`
to diff the non-tick columns and the whole LPT assignment across two runs.

The tables above and the figures come from the same run directory:

```
python benchmarks/spec_matching/plot_load_balancer.py results --table
```

It writes `load_balancer_latency.txt` (the text table this file is built from) plus three families
of figure: `lbp_d*_p*_T*_k*.png`, one histogram of `fallback` against `sparse_k` per cell;
`lbp_ksweep_d*_p*_T*.png`, the same pair as small multiples with one panel per `k` on shared bins;
and `lbp_modules_d*_p*_T*.png`, the six timed regions as horizontal bars, one panel per `k`. Pass
`--subtract-overhead` for the per-column timer correction described above, `--log-y` to bring the
escalation tail up, and `--d` / `--T` / `--k` to cut the sweep down.

For the escalation rate `q` and the structure of `H` on their own — no timers, no decoder output,
a Wilson interval on `q` and a per-component `COMPLETE`/`TRUNCATED` breakdown — use the companion
profiler instead:

```
cmake --build build --target sparse_graph_stats -j4
build/sparse_graph_stats --d 17,21,25,31 --p 0.0005,0.001,0.003 --T 1.5,2 \
    --shots 1000000 --seed 20260907 --out results/sparse_graph_stats
```

The run directory this file was written from is not in the repository — `shots.csv` alone is 8.2 GB.
It lives outside the tree, alongside it, in `../results`.
