# Latency speedups

> **What this file is now.** Every number below was produced by `load_balancer_profiler` at git
> `0155f2c`, on Apple silicon under macOS, and measures the **k-core critical path of the decomposed
> `H` solve** against stock exact decode on `G`. It replaces the resolver-era sweep this file used to
> carry. In those rows `k` was a component-size cut for a lookup-table resolver that no longer
> exists, and every row naming a `k` other than 0 decoded a different decoder output; the file said
> so at the top and the rows are not reproduced here. They are in git at `dccfcee`, which is also
> where the reason they are wrong is written down. **In this file `k` is the number of solver
> cores.** Nothing else about the two sweeps is comparable either — different machine, different
> timer, different baseline series, no contamination filter.

Per-shot latency of the spec-matching decoder when the solve on the sparsified graph `H` is
decomposed across `k` cores, against stock exact sparse blossom on `G`, over a 48-point sweep of a
rotated surface-code memory-X experiment.

## What was measured

The profiler logs three latencies per shot. `run.log` carries the definitions; they are repeated
here because every column in this file is one of them:

```
sparse_k = uf + balance + max_c(build_c + solve_c + extract_c) + combine
fallback = stock exact decode on G, from the full syndrome, on its own instance, run on EVERY shot
system   = escalated ? fallback : min(fallback, sparse_k)
```

* **`uf`** — union-find over `H`'s defect-defect edges: component id per defect, component sizes.
* **`balance`** — cost lookup, LPT assignment, and the scatter of defects, edges and boundary edges
  into the per-core slices. The balancer's own statistics are read after the region closes.
* **`build_c`** — sub-`H_c` construction from core `c`'s slice, `BallMwpm::rebuild`, and the
  detection-event fill.
* **`solve_c`** — `process_timeline_until_horizon` on that core's sub-graph, and nothing else.
* **`extract_c`** — `extract_only_to_obs` alone, on a `COMPLETE` core. `reset_for_next_shot` still
  runs, immediately after, but outside the region: this shot's answer does not depend on it, so a
  double-buffered decoder would not pay it here.
* **`combine`** — XOR of the `k` partial observable masks, sum of the partial weights, OR of the
  truncation flags.

The three build-solve-extract regions are charged as a **max over cores, never a sum**; `uf`,
`balance` and `combine` are serial and are charged in full. Detection-event generation and the ball
graph build are the *input* and sit outside every timed region.

A shot **escalates** when any core's timeline is `TRUNCATED`, which is the same as any component
truncating, because a core's timeline is the disjoint union of its components'. The escalation counts
are identical at `k = 5` and `k = 9` in all 24 `(d, p, T)` corpora — the decomposition does not move
the escalating set, which is what the independence property predicts.

`speedup` is `mean(fallback) / mean(system)`, over every shot including the escalating ones.

### The model, and what it does not charge

`k_core_critical_path_model`: **no thread is ever spawned.** The `k` loads are solved one after
another on the profiling thread, each on its own pre-allocated `Mwpm` instance, and the per-core
figure is the max. Four consequences travel with every number here:

1. **No cost is modelled for moving a core's edge slice to it.** The scatter happens once, serially,
   inside `balance`, on the profiling thread. Conservative in the sense that a real fan-out would add
   to `balance`, not to the per-core regions.
2. **A load solved straight after the previous one sees a warmer cache than a real core would.** This
   biases `sparse_k` low and the speedup high, by an amount this experiment cannot bound.
3. **The fallback core and the `k` solver cores are assumed to start at the same instant**, when the
   edge list is ready, and to share nothing.
4. **Rows are raw, uncorrected tick deltas.** `sparse_k` sums six bracketed intervals against
   `fallback`'s one, so it carries about five extra timer reads — 155.117 ns each, measured over 10^6
   back-to-back calls. `plot_load_balancer.py --subtract-overhead` does the per-column correction;
   the tables here are raw, so `sparse_k` is charged ~0.8 µs it would not pay in a real decoder.

### Sweep

| | |
|---|---|
| circuit | `surface_code:rotated_memory_x`, `rounds = d`, `p` on Clifford / reset / measurement |
| distances `d` | 17, 21, 25, 31 |
| error rates `p` | 0.0005, 0.001, 0.003 |
| horizon `T` | 1.5, 2.0 (units of the median discretised edge of `G`) |
| cores `k` | 5, 9 |
| cost model | `cost(s) = s^alpha`, `alpha = 1` — component **size** only, tabulated once per `(d, p)` |
| balancer | LPT: components in descending size, each to the least-loaded core |
| shots | 100,000 per cell, 48 cells, 4,800,000 rows; 1,000 warm-up per cell |
| seed | 20260906 — shot `i` of a `(d, p)` corpus is the same shot across the `T` and `k` sweeps |
| timer | `thread_cputime_macos`, thread-scoped, 155.117 ns/read |
| build | `-O3`, stim v1.16.0, macOS |

Warm-up shots run the identical workload — same builds, solves, extractions and fallback, in the
same order, on the same buffers — and **read no clock at all**: `run_shot` is templated on the
switch, so the reads compile out of the warm-up instantiation. There is no contamination probe and
no shot is dropped; all 100,000 rows per cell are in every number below.

### One caveat that is not the decoder's

`fallback` is the internal control: the same shots, decoded the same way, in every cell of a
`(d, p, T)` corpus. It should not move with `k`. **It moves.** The cells ran back to back over
70 minutes on a laptop, and at `d = 21` the control drifts up to 36% between the `k = 5` and `k = 9`
cells:

| d | p | T | `fallback` k9/k5 | `sparse_k` k9/k5 | speedup k=5 → k=9 |
|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 1.020 | 0.976 | 2.14 → 2.23 |
| 17 | 0.0005 | 2 | 1.012 | 1.018 | 1.67 → 1.66 |
| 17 | 0.001 | 1.5 | 1.135 | 0.970 | 1.93 → 2.24 |
| 17 | 0.001 | 2 | 1.101 | 1.025 | 1.28 → 1.35 |
| 17 | 0.003 | 1.5 | 1.091 | 1.042 | 1.11 → 1.13 |
| 17 | 0.003 | 2 | 1.010 | 1.008 | 1.00 → 1.00 |
| 21 | 0.0005 | 1.5 | 1.363 | 1.011 | 2.61 → 3.50 |
| 21 | 0.0005 | 2 | 1.299 | 1.043 | 2.26 → 2.81 |
| 21 | 0.001 | 1.5 | 1.253 | 0.936 | 2.67 → 3.50 |
| 21 | 0.001 | 2 | 1.219 | 1.062 | 1.66 → 1.89 |
| 21 | 0.003 | 1.5 | 1.154 | 1.047 | 1.11 → 1.16 |
| 21 | 0.003 | 2 | 1.006 | 1.006 | 1.00 → 1.00 |
| 25 | 0.0005 | 1.5 | 1.186 | 0.940 | 3.89 → 4.86 |
| 25 | 0.0005 | 2 | 1.169 | 0.948 | 3.31 → 4.08 |
| 25 | 0.001 | 1.5 | 1.142 | 0.867 | 3.64 → 4.62 |
| 25 | 0.001 | 2 | 1.101 | 1.035 | 2.36 → 2.51 |
| 25 | 0.003 | 1.5 | 1.067 | 1.021 | 1.12 → 1.15 |
| 25 | 0.003 | 2 | 1.007 | 1.005 | 1.00 → 1.00 |
| 31 | 0.0005 | 1.5 | 1.024 | 0.881 | 4.81 → 5.53 |
| 31 | 0.0005 | 2 | 1.028 | 0.858 | 4.27 → 5.12 |
| 31 | 0.001 | 1.5 | 1.025 | 0.804 | 3.98 → 4.80 |
| 31 | 0.001 | 2 | 1.055 | 1.000 | 2.71 → 2.86 |
| 31 | 0.003 | 1.5 | 1.013 | 1.010 | 1.09 → 1.09 |
| 31 | 0.003 | 2 | 0.995 | 1.000 | 1.00 → 1.00 |

So **read the `k = 5` → `k = 9` speedup gain at `d = 21` and `d = 25` as mostly drift, not as
balancer.** The `sparse_k` column is the honest one for the `k` question, and it is the one to
believe at `d = 31`, where the control is stable to 2–5% and `sparse_k` genuinely drops 12–20%.
Every absolute latency in this file is a laptop number and none of them should be quoted as a
hardware figure; the ratios within a cell are what the sweep supports.

## Results

Speedup ranges from **1.00x** (every `p = 0.003, T = 2` cell — see below) to **5.53x**
(`d = 31, p = 0.0005, T = 1.5, k = 9`). The pattern:

* **It grows with `d` and collapses with `p`,** as the previous sweep found and for the same reason:
  the fallback scales with all of `G`, the critical path scales with the largest *component* of `H`,
  and at `p = 0.003` `H` stops being fragmented. At `p = 0.0005, T = 1.5, k = 9` the speedup goes
  2.23 → 3.50 → 4.86 → 5.53 across `d = 17, 21, 25, 31`; at `p = 0.003, T = 2` it is 1.00 at every
  distance.
* **At `p = 0.003, T = 2` the decomposition is a strict loss and `min` is the only thing saving it.**
  `H` is one component holding almost every defect — 1382.7 of 1391.6 at `d = 31` — so there is
  nothing to spread, `sparse_k` runs at 1.7x the fallback (1016.7 µs against 613.4 µs), it beats the
  fallback on **0.0%** of shots, and `system` is `fallback` on essentially every shot. The 1.00 is
  not a wash between two comparable paths; it is the concurrent structure refusing to lose.
* **The bottleneck at high `p` is `build`, not `solve`.** `build_crit` is the largest module in every
  cell above `p = 0.0005`, and at `p = 0.003, T = 2, d = 31` it is 613.6 µs against `solve_crit`'s
  310.3 µs — sub-`H_c` construction and `BallMwpm::rebuild` cost twice what blossom costs on the same
  graph, and alone they exceed the whole stock decode. Whatever `k` buys, it does not touch this.
* **The serial preamble is the ceiling on `k`.** `uf + balance + combine` do not parallelize, and at
  low `p` they are most of the critical path: 51% of `sparse_k` at `d = 31, p = 0.0005, T = 1.5,
  k = 9`, and 38–47% across the other low-`p`, `k = 9` cells at `d = 25` and `d = 31`. Going from 5 cores to 9
  moves work out of `max_c` and **into** `balance` — at `d = 31, p = 0.0005, T = 1.5`, `balance`
  goes 3.96 → 5.42 µs while the per-core max goes 12.55 → 8.56 µs. That is why the `sparse_k` gain
  from doubling the cores is ~12–20% and not ~2x, and it is a real ceiling rather than a
  measurement artefact.
* **LPT balances well exactly when there is something to balance.** Mean imbalance
  (`k * pred_max / pred_sum`; 1.0 is perfect) is 1.00–1.06 wherever components outnumber cores by a
  wide margin, and degrades to ≈ `k` — 8.94 at `k = 9` — when `H` is a single giant component and one
  core holds everything. The balancer is not the problem at `p = 0.003`; the absence of components
  is.
* **`T = 1.5` beats `T = 2` on the mean everywhere, and pays in escalation**, from 0.01% to 1.2% at
  `p = 0.0005` and from 2.6% to 28.9% at `p = 0.003, d = 31`. Because an escalating shot costs a lost
  race and not a second decode, the tighter horizon still wins at every point.

### Speedup table

All latencies in microseconds. `escalated` is the share of shots on which some core truncated;
`mean x` is `mean(fallback) / mean(system)`, the percentile columns are ratios of the corresponding
percentiles. Read the `k = 9` rows against the drift table above before attributing a gain to `k`.

| d | p | T | k | escalated | fallback mean | sparse_k mean | system mean | **mean x** | p50 x | p99 x | p99.9 x |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 0.321% | 9.05 | 4.20 | 4.24 | **2.14** | 2.19 | 2.00 | 1.21 |
| 17 | 0.0005 | 1.5 | 9 | 0.321% | 9.23 | 4.10 | 4.13 | **2.23** | 2.30 | 2.01 | 1.21 |
| 17 | 0.0005 | 2 | 5 | 0.002% | 8.85 | 5.35 | 5.30 | **1.67** | 1.77 | 1.39 | 1.36 |
| 17 | 0.0005 | 2 | 9 | 0.002% | 8.95 | 5.44 | 5.39 | **1.66** | 1.75 | 1.39 | 1.38 |
| 17 | 0.001 | 1.5 | 5 | 1.382% | 15.33 | 7.75 | 7.94 | **1.93** | 2.01 | 1.34 | 1.14 |
| 17 | 0.001 | 1.5 | 9 | 1.382% | 17.40 | 7.52 | 7.75 | **2.24** | 2.42 | 1.36 | 1.13 |
| 17 | 0.001 | 2 | 5 | 0.021% | 16.43 | 13.96 | 12.82 | **1.28** | 1.32 | 1.10 | 1.14 |
| 17 | 0.001 | 2 | 9 | 0.021% | 18.09 | 14.30 | 13.42 | **1.35** | 1.40 | 1.11 | 1.18 |
| 17 | 0.003 | 1.5 | 5 | 7.721% | 61.45 | 61.73 | 55.12 | **1.11** | 1.08 | 1.01 | 1.01 |
| 17 | 0.003 | 1.5 | 9 | 7.721% | 67.05 | 64.31 | 59.09 | **1.13** | 1.10 | 1.02 | 1.02 |
| 17 | 0.003 | 2 | 5 | 0.627% | 64.79 | 128.36 | 64.78 | **1.00** | 1.00 | 1.00 | 1.00 |
| 17 | 0.003 | 2 | 9 | 0.627% | 65.46 | 129.38 | 65.46 | **1.00** | 1.00 | 1.00 | 1.00 |
| 21 | 0.0005 | 1.5 | 5 | 0.511% | 17.12 | 6.48 | 6.56 | **2.61** | 2.66 | 2.38 | 1.23 |
| 21 | 0.0005 | 1.5 | 9 | 0.511% | 23.33 | 6.55 | 6.67 | **3.50** | 3.63 | 2.99 | 1.20 |
| 21 | 0.0005 | 2 | 5 | 0.010% | 18.09 | 8.01 | 8.00 | **2.26** | 2.35 | 1.82 | 1.66 |
| 21 | 0.0005 | 2 | 9 | 0.010% | 23.50 | 8.35 | 8.35 | **2.81** | 2.98 | 2.08 | 1.82 |
| 21 | 0.001 | 1.5 | 5 | 2.053% | 37.11 | 13.26 | 13.91 | **2.67** | 2.82 | 1.23 | 1.14 |
| 21 | 0.001 | 1.5 | 9 | 2.053% | 46.51 | 12.41 | 13.28 | **3.50** | 3.89 | 1.22 | 1.11 |
| 21 | 0.001 | 2 | 5 | 0.038% | 38.61 | 23.54 | 23.20 | **1.66** | 1.84 | 1.11 | 1.12 |
| 21 | 0.001 | 2 | 9 | 0.038% | 47.06 | 25.00 | 24.89 | **1.89** | 2.09 | 1.16 | 1.13 |
| 21 | 0.003 | 1.5 | 5 | 12.632% | 141.57 | 134.28 | 127.19 | **1.11** | 1.06 | 1.01 | 1.01 |
| 21 | 0.003 | 1.5 | 9 | 12.632% | 163.43 | 140.55 | 140.73 | **1.16** | 1.10 | 1.01 | 1.01 |
| 21 | 0.003 | 2 | 5 | 0.992% | 144.58 | 272.03 | 144.58 | **1.00** | 1.00 | 1.00 | 1.00 |
| 21 | 0.003 | 2 | 9 | 0.992% | 145.47 | 273.74 | 145.47 | **1.00** | 1.00 | 1.00 | 1.00 |
| 25 | 0.0005 | 1.5 | 5 | 0.747% | 42.48 | 10.63 | 10.93 | **3.89** | 4.03 | 3.19 | 1.15 |
| 25 | 0.0005 | 1.5 | 9 | 0.747% | 50.38 | 9.99 | 10.36 | **4.86** | 5.11 | 3.91 | 1.25 |
| 25 | 0.0005 | 2 | 5 | 0.012% | 42.17 | 12.74 | 12.75 | **3.31** | 3.38 | 2.74 | 1.97 |
| 25 | 0.0005 | 2 | 9 | 0.012% | 49.30 | 12.08 | 12.08 | **4.08** | 4.29 | 3.01 | 2.31 |
| 25 | 0.001 | 1.5 | 5 | 3.068% | 87.62 | 21.74 | 24.05 | **3.64** | 4.05 | 1.13 | 1.07 |
| 25 | 0.001 | 1.5 | 9 | 3.068% | 100.02 | 18.86 | 21.64 | **4.62** | 5.41 | 1.13 | 1.08 |
| 25 | 0.001 | 2 | 5 | 0.062% | 88.76 | 37.51 | 37.54 | **2.36** | 2.67 | 1.31 | 1.17 |
| 25 | 0.001 | 2 | 9 | 0.062% | 97.68 | 38.81 | 38.86 | **2.51** | 2.78 | 1.37 | 1.19 |
| 25 | 0.003 | 1.5 | 5 | 18.468% | 291.19 | 254.29 | 259.15 | **1.12** | 1.06 | 1.01 | 1.00 |
| 25 | 0.003 | 1.5 | 9 | 18.468% | 310.59 | 259.52 | 269.65 | **1.15** | 1.09 | 1.01 | 1.01 |
| 25 | 0.003 | 2 | 5 | 1.478% | 284.62 | 493.57 | 284.62 | **1.00** | 1.00 | 1.00 | 1.00 |
| 25 | 0.003 | 2 | 9 | 1.478% | 286.56 | 496.15 | 286.56 | **1.00** | 1.00 | 1.00 | 1.00 |
| 31 | 0.0005 | 1.5 | 5 | 1.197% | 99.73 | 19.67 | 20.75 | **4.81** | 5.12 | 1.28 | 1.11 |
| 31 | 0.0005 | 1.5 | 9 | 1.197% | 102.12 | 17.33 | 18.47 | **5.53** | 5.95 | 1.28 | 1.11 |
| 31 | 0.0005 | 2 | 5 | 0.010% | 98.76 | 23.12 | 23.13 | **4.27** | 4.32 | 3.91 | 2.44 |
| 31 | 0.0005 | 2 | 9 | 0.010% | 101.54 | 19.84 | 19.85 | **5.12** | 5.24 | 4.09 | 2.65 |
| 31 | 0.001 | 1.5 | 5 | 4.893% | 197.17 | 41.40 | 49.54 | **3.98** | 4.79 | 1.07 | 1.04 |
| 31 | 0.001 | 1.5 | 9 | 4.893% | 202.15 | 33.30 | 42.09 | **4.80** | 6.12 | 1.07 | 1.05 |
| 31 | 0.001 | 2 | 5 | 0.101% | 189.55 | 69.68 | 69.83 | **2.71** | 3.06 | 1.51 | 1.20 |
| 31 | 0.001 | 2 | 9 | 0.101% | 200.05 | 69.70 | 69.86 | **2.86** | 3.13 | 1.56 | 1.24 |
| 31 | 0.003 | 1.5 | 5 | 28.903% | 624.14 | 544.26 | 572.18 | **1.09** | 1.05 | 1.01 | 1.00 |
| 31 | 0.003 | 1.5 | 9 | 28.903% | 632.33 | 549.96 | 579.10 | **1.09** | 1.05 | 1.01 | 1.00 |
| 31 | 0.003 | 2 | 5 | 2.620% | 613.35 | 1016.70 | 613.35 | **1.00** | 1.00 | 1.00 | 1.00 |
| 31 | 0.003 | 2 | 9 | 2.620% | 610.43 | 1016.88 | 610.43 | **1.00** | 1.00 | 1.00 | 1.00 |


### Tails, in microseconds

The mean is not what sizes a syndrome buffer — blocks arrive on the hardware's clock whether or not
the last one is decoded — and the tail is where the two series converge. The convergence is the
expected shape, not a defect: a shot in the tail is a shot with a large component in `H`, which is
exactly the shot that truncates, and on a truncating shot `system` *is* `fallback` by construction.
The tail holds up only where escalation is rare enough that the 99.9th percentile is still a
non-escalating shot — `d = 31, p = 0.0005, T = 2, k = 9` keeps 2.65x at p99.9, against 1.00x at
`p = 0.003`.

The `system` max is worth reading beside the `fallback` max: at `d = 17, p = 0.0005, T = 1.5, k = 5`
the worst `fallback` shot is 56.5 µs and the worst `system` shot 30.1 µs, because `system` takes the
min and the fallback's worst shot is not the sparse path's worst shot.

| d | p | T | k | fb p50 | sys p50 | fb p99 | sys p99 | fb p99.9 | sys p99.9 | fb max | sys max |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 8.75 | 4.00 | 17.25 | 8.62 | 21.79 | 18.08 | 56.5 | 30.1 |
| 17 | 0.0005 | 1.5 | 9 | 8.92 | 3.88 | 17.58 | 8.75 | 22.17 | 18.29 | 67.6 | 33.7 |
| 17 | 0.0005 | 2 | 5 | 8.54 | 4.83 | 16.88 | 12.17 | 22.08 | 16.21 | 51.2 | 32.8 |
| 17 | 0.0005 | 2 | 9 | 8.62 | 4.92 | 17.17 | 12.37 | 22.58 | 16.33 | 57.2 | 33.3 |
| 17 | 0.001 | 1.5 | 5 | 14.83 | 7.37 | 27.29 | 20.42 | 35.67 | 31.17 | 60.8 | 58.1 |
| 17 | 0.001 | 1.5 | 9 | 16.88 | 6.96 | 30.38 | 22.38 | 38.38 | 34.00 | 66.0 | 60.7 |
| 17 | 0.001 | 2 | 5 | 15.92 | 12.04 | 28.79 | 26.29 | 37.25 | 32.75 | 63.9 | 63.3 |
| 17 | 0.001 | 2 | 9 | 17.54 | 12.50 | 31.17 | 28.17 | 40.88 | 34.67 | 69.5 | 67.0 |
| 17 | 0.003 | 1.5 | 5 | 59.88 | 55.21 | 99.33 | 97.92 | 123.54 | 122.50 | 352.5 | 352.5 |
| 17 | 0.003 | 1.5 | 9 | 65.50 | 59.59 | 106.42 | 104.46 | 132.67 | 130.00 | 367.8 | 367.8 |
| 17 | 0.003 | 2 | 5 | 63.25 | 63.25 | 103.62 | 103.62 | 129.04 | 129.04 | 382.2 | 382.2 |
| 17 | 0.003 | 2 | 9 | 63.92 | 63.92 | 104.79 | 104.79 | 129.46 | 129.46 | 362.8 | 362.8 |
| 21 | 0.0005 | 1.5 | 5 | 16.75 | 6.29 | 28.17 | 11.83 | 35.88 | 29.17 | 62.7 | 45.5 |
| 21 | 0.0005 | 1.5 | 9 | 23.00 | 6.33 | 36.50 | 12.21 | 45.38 | 37.67 | 75.3 | 58.5 |
| 21 | 0.0005 | 2 | 5 | 17.71 | 7.54 | 29.62 | 16.29 | 37.83 | 22.83 | 65.2 | 47.6 |
| 21 | 0.0005 | 2 | 9 | 23.12 | 7.75 | 36.71 | 17.67 | 44.54 | 24.46 | 73.4 | 55.2 |
| 21 | 0.001 | 1.5 | 5 | 36.58 | 12.96 | 55.42 | 45.12 | 68.54 | 60.33 | 97.4 | 90.7 |
| 21 | 0.001 | 1.5 | 9 | 46.00 | 11.83 | 66.96 | 55.08 | 79.88 | 71.79 | 112.2 | 112.2 |
| 21 | 0.001 | 2 | 5 | 38.08 | 20.67 | 57.38 | 51.46 | 68.12 | 61.04 | 95.2 | 91.5 |
| 21 | 0.001 | 2 | 9 | 46.50 | 22.29 | 67.92 | 58.75 | 79.75 | 70.62 | 119.3 | 103.7 |
| 21 | 0.003 | 1.5 | 5 | 139.83 | 131.67 | 195.46 | 194.04 | 226.17 | 223.92 | 314.4 | 314.4 |
| 21 | 0.003 | 1.5 | 9 | 161.71 | 146.83 | 219.88 | 217.38 | 251.00 | 248.54 | 339.6 | 339.6 |
| 21 | 0.003 | 2 | 5 | 142.83 | 142.83 | 198.62 | 198.62 | 229.12 | 229.12 | 312.8 | 312.8 |
| 21 | 0.003 | 2 | 9 | 143.67 | 143.67 | 199.71 | 199.71 | 228.96 | 228.96 | 318.5 | 318.5 |
| 25 | 0.0005 | 1.5 | 5 | 42.12 | 10.46 | 60.50 | 18.96 | 70.38 | 61.17 | 94.9 | 93.2 |
| 25 | 0.0005 | 1.5 | 9 | 50.00 | 9.79 | 70.58 | 18.04 | 86.33 | 69.33 | 114.6 | 92.6 |
| 25 | 0.0005 | 2 | 5 | 41.75 | 12.33 | 60.71 | 22.17 | 72.29 | 36.71 | 141.6 | 75.8 |
| 25 | 0.0005 | 2 | 9 | 48.96 | 11.42 | 68.83 | 22.87 | 79.71 | 34.50 | 105.7 | 85.8 |
| 25 | 0.001 | 1.5 | 5 | 87.00 | 21.46 | 115.75 | 102.46 | 129.79 | 121.00 | 162.1 | 162.1 |
| 25 | 0.001 | 1.5 | 9 | 99.42 | 18.38 | 129.71 | 115.21 | 145.08 | 134.46 | 177.3 | 177.3 |
| 25 | 0.001 | 2 | 5 | 88.21 | 33.04 | 117.42 | 89.58 | 131.67 | 112.88 | 162.1 | 161.9 |
| 25 | 0.001 | 2 | 9 | 97.12 | 34.96 | 127.25 | 93.04 | 142.29 | 119.96 | 194.8 | 172.5 |
| 25 | 0.003 | 1.5 | 5 | 289.21 | 271.75 | 365.88 | 362.88 | 403.46 | 401.58 | 550.7 | 550.7 |
| 25 | 0.003 | 1.5 | 9 | 308.54 | 282.12 | 388.08 | 383.42 | 427.17 | 423.58 | 580.3 | 580.3 |
| 25 | 0.003 | 2 | 5 | 282.42 | 282.42 | 361.08 | 361.08 | 403.50 | 403.50 | 538.5 | 538.5 |
| 25 | 0.003 | 2 | 9 | 284.46 | 284.46 | 361.29 | 361.29 | 401.38 | 401.38 | 541.8 | 541.8 |
| 31 | 0.0005 | 1.5 | 5 | 99.33 | 19.42 | 127.46 | 99.21 | 141.38 | 127.25 | 176.5 | 157.7 |
| 31 | 0.0005 | 1.5 | 9 | 101.71 | 17.08 | 130.17 | 101.54 | 143.33 | 129.58 | 182.1 | 165.2 |
| 31 | 0.0005 | 2 | 5 | 98.38 | 22.79 | 126.62 | 32.42 | 142.71 | 58.37 | 168.2 | 135.5 |
| 31 | 0.0005 | 2 | 9 | 101.12 | 19.29 | 129.58 | 31.67 | 144.50 | 54.54 | 184.7 | 139.3 |
| 31 | 0.001 | 1.5 | 5 | 196.62 | 41.04 | 238.08 | 222.08 | 256.42 | 246.17 | 332.0 | 332.0 |
| 31 | 0.001 | 1.5 | 9 | 201.54 | 32.92 | 244.04 | 227.12 | 263.67 | 251.25 | 315.3 | 288.1 |
| 31 | 0.001 | 2 | 5 | 188.88 | 61.67 | 231.42 | 153.33 | 252.42 | 209.92 | 311.0 | 271.9 |
| 31 | 0.001 | 2 | 9 | 199.25 | 63.59 | 245.21 | 157.50 | 267.96 | 216.17 | 360.9 | 284.4 |
| 31 | 0.003 | 1.5 | 5 | 621.62 | 594.25 | 734.54 | 730.08 | 787.42 | 785.25 | 1161.0 | 1161.0 |
| 31 | 0.003 | 1.5 | 9 | 629.92 | 601.33 | 744.38 | 740.21 | 797.04 | 793.71 | 1172.1 | 1172.1 |
| 31 | 0.003 | 2 | 5 | 611.08 | 611.08 | 721.38 | 721.38 | 773.92 | 773.92 | 1132.6 | 1132.6 |
| 31 | 0.003 | 2 | 9 | 608.21 | 608.21 | 716.67 | 716.67 | 768.00 | 768.00 | 1180.9 | 1180.9 |


### Where the critical path goes

Module means in microseconds, from `agg.json`'s 128-bit accumulators. The six sum to the `sparse_k`
mean exactly, by the definition at the top. `serial share` is `(uf + balance + combine) / sparse_k` —
the part of the critical path that does not shrink with `k`, no matter how many cores are added.

| d | p | T | k | uf | balance | build | solve | extract | combine | sparse_k | serial share |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 0.54 | 0.68 | 1.17 | 1.31 | 0.33 | 0.163 | 4.20 | 33% |
| 17 | 0.0005 | 1.5 | 9 | 0.54 | 0.92 | 1.06 | 1.11 | 0.30 | 0.175 | 4.10 | 40% |
| 17 | 0.0005 | 2 | 5 | 0.62 | 0.64 | 1.87 | 1.69 | 0.37 | 0.162 | 5.35 | 27% |
| 17 | 0.0005 | 2 | 9 | 0.62 | 0.83 | 1.83 | 1.61 | 0.37 | 0.177 | 5.44 | 30% |
| 17 | 0.001 | 1.5 | 5 | 1.05 | 1.06 | 2.33 | 2.66 | 0.50 | 0.162 | 7.75 | 29% |
| 17 | 0.001 | 1.5 | 9 | 1.06 | 1.39 | 2.18 | 2.28 | 0.44 | 0.176 | 7.52 | 35% |
| 17 | 0.001 | 2 | 5 | 1.25 | 1.02 | 5.97 | 4.80 | 0.75 | 0.162 | 13.96 | 17% |
| 17 | 0.001 | 2 | 9 | 1.27 | 1.22 | 6.08 | 4.80 | 0.75 | 0.175 | 14.30 | 19% |
| 17 | 0.003 | 1.5 | 5 | 3.66 | 2.54 | 26.82 | 25.59 | 2.96 | 0.163 | 61.73 | 10% |
| 17 | 0.003 | 1.5 | 9 | 3.72 | 2.75 | 28.67 | 25.99 | 2.98 | 0.205 | 64.31 | 10% |
| 17 | 0.003 | 2 | 5 | 4.75 | 3.69 | 69.08 | 45.93 | 4.76 | 0.163 | 128.36 | 7% |
| 17 | 0.003 | 2 | 9 | 4.75 | 3.73 | 69.99 | 45.98 | 4.75 | 0.168 | 129.38 | 7% |
| 21 | 0.0005 | 1.5 | 5 | 0.90 | 1.18 | 1.71 | 2.07 | 0.46 | 0.162 | 6.48 | 35% |
| 21 | 0.0005 | 1.5 | 9 | 0.97 | 1.65 | 1.69 | 1.64 | 0.36 | 0.236 | 6.55 | 44% |
| 21 | 0.0005 | 2 | 5 | 1.06 | 1.08 | 2.73 | 2.49 | 0.48 | 0.162 | 8.01 | 29% |
| 21 | 0.0005 | 2 | 9 | 1.12 | 1.44 | 2.88 | 2.24 | 0.43 | 0.239 | 8.35 | 34% |
| 21 | 0.001 | 1.5 | 5 | 1.95 | 1.87 | 4.10 | 4.43 | 0.75 | 0.167 | 13.26 | 30% |
| 21 | 0.001 | 1.5 | 9 | 2.05 | 2.59 | 3.62 | 3.33 | 0.56 | 0.241 | 12.41 | 39% |
| 21 | 0.001 | 2 | 5 | 2.31 | 1.77 | 10.56 | 7.62 | 1.11 | 0.183 | 23.54 | 18% |
| 21 | 0.001 | 2 | 9 | 2.42 | 2.19 | 11.44 | 7.61 | 1.09 | 0.244 | 25.00 | 19% |
| 21 | 0.003 | 1.5 | 5 | 7.14 | 4.82 | 63.12 | 53.19 | 5.84 | 0.169 | 134.28 | 9% |
| 21 | 0.003 | 1.5 | 9 | 7.23 | 5.22 | 68.09 | 53.90 | 5.88 | 0.246 | 140.55 | 9% |
| 21 | 0.003 | 2 | 5 | 9.35 | 7.10 | 155.63 | 90.42 | 9.36 | 0.168 | 272.03 | 6% |
| 21 | 0.003 | 2 | 9 | 9.36 | 7.15 | 157.26 | 90.42 | 9.38 | 0.174 | 273.74 | 6% |
| 25 | 0.0005 | 1.5 | 5 | 1.52 | 1.95 | 2.97 | 3.35 | 0.64 | 0.190 | 10.63 | 35% |
| 25 | 0.0005 | 1.5 | 9 | 1.62 | 2.81 | 2.49 | 2.35 | 0.46 | 0.260 | 9.99 | 47% |
| 25 | 0.0005 | 2 | 5 | 1.78 | 1.79 | 4.45 | 3.88 | 0.65 | 0.195 | 12.74 | 30% |
| 25 | 0.0005 | 2 | 9 | 1.88 | 2.43 | 4.02 | 2.98 | 0.51 | 0.261 | 12.08 | 38% |
| 25 | 0.001 | 1.5 | 5 | 3.31 | 3.17 | 6.84 | 7.04 | 1.15 | 0.230 | 21.74 | 31% |
| 25 | 0.001 | 1.5 | 9 | 3.48 | 4.34 | 5.19 | 4.84 | 0.76 | 0.250 | 18.86 | 43% |
| 25 | 0.001 | 2 | 5 | 3.91 | 2.94 | 17.17 | 11.62 | 1.63 | 0.231 | 37.51 | 19% |
| 25 | 0.001 | 2 | 9 | 4.10 | 3.61 | 17.96 | 11.32 | 1.57 | 0.253 | 38.81 | 21% |
| 25 | 0.003 | 1.5 | 5 | 12.48 | 8.11 | 124.44 | 98.56 | 10.48 | 0.231 | 254.29 | 8% |
| 25 | 0.003 | 1.5 | 9 | 12.51 | 8.89 | 127.98 | 99.28 | 10.61 | 0.255 | 259.52 | 8% |
| 25 | 0.003 | 2 | 5 | 16.55 | 12.37 | 290.22 | 157.41 | 16.83 | 0.199 | 493.57 | 6% |
| 25 | 0.003 | 2 | 9 | 16.53 | 12.45 | 292.85 | 157.30 | 16.81 | 0.204 | 496.15 | 6% |
| 31 | 0.0005 | 1.5 | 5 | 2.94 | 3.96 | 5.64 | 5.85 | 1.06 | 0.229 | 19.67 | 36% |
| 31 | 0.0005 | 1.5 | 9 | 3.04 | 5.42 | 3.96 | 3.92 | 0.68 | 0.302 | 17.33 | 51% |
| 31 | 0.0005 | 2 | 5 | 3.39 | 3.62 | 8.21 | 6.60 | 1.07 | 0.232 | 23.12 | 31% |
| 31 | 0.0005 | 2 | 9 | 3.49 | 4.70 | 6.13 | 4.51 | 0.70 | 0.302 | 19.84 | 43% |
| 31 | 0.001 | 1.5 | 5 | 6.41 | 6.23 | 13.72 | 12.71 | 2.08 | 0.246 | 41.40 | 31% |
| 31 | 0.001 | 1.5 | 9 | 6.47 | 8.19 | 8.98 | 8.13 | 1.27 | 0.261 | 33.30 | 45% |
| 31 | 0.001 | 2 | 5 | 7.74 | 5.76 | 32.84 | 20.31 | 2.80 | 0.248 | 69.68 | 20% |
| 31 | 0.001 | 2 | 9 | 7.79 | 6.84 | 32.78 | 19.38 | 2.65 | 0.263 | 69.70 | 21% |
| 31 | 0.003 | 1.5 | 5 | 24.15 | 15.76 | 274.08 | 207.81 | 22.22 | 0.244 | 544.26 | 7% |
| 31 | 0.003 | 1.5 | 9 | 24.24 | 16.85 | 277.49 | 208.78 | 22.32 | 0.274 | 549.96 | 8% |
| 31 | 0.003 | 2 | 5 | 32.37 | 24.23 | 613.55 | 310.25 | 36.08 | 0.227 | 1016.70 | 6% |
| 31 | 0.003 | 2 | 9 | 32.40 | 24.36 | 614.15 | 309.67 | 36.08 | 0.233 | 1016.88 | 6% |

### Structure of H, and how well it balances

Means over the same shots, all read outside every timed region. `crit size` and `crit comps` describe
the core that turned out to be critical — the argmax of measured ticks, so unlike the rest of this
table they are not in the deterministic set. `imbalance` is `k * pred_max / pred_sum` from the LPT
prediction, where 1.0 is perfect and `k` is one core holding everything. `sparse_k wins` is the share
of shots on which the decomposed path finished before the fallback — the share on which the `min` in
`system` selects the sparse side.

This table is the explanation for the ones above, not another measurement of them. The three regimes
are visible in `components` against `largest`: many small components (`d = 31, p = 0.0005`: 89.8
components, largest 10.8) where the decomposition wins outright; a transition (`p = 0.001, T = 2`);
and one giant component (`p = 0.003, T = 2`: 4.5 components, largest 1382.7) where there is nothing
to decompose.

| d | p | T | k | defects | components | largest | crit size | crit comps | imbalance | sparse_k wins |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 5 | 38.3 | 14.92 | 6.17 | 8.47 | 2.37 | 1.14 | 99.8% |
| 17 | 0.0005 | 1.5 | 9 | 38.3 | 14.92 | 6.17 | 6.33 | 1.20 | 1.54 | 99.8% |
| 17 | 0.0005 | 2 | 5 | 38.3 | 11.32 | 9.69 | 10.34 | 1.47 | 1.36 | 96.9% |
| 17 | 0.0005 | 2 | 9 | 38.3 | 11.32 | 9.69 | 9.60 | 1.02 | 2.25 | 96.2% |
| 17 | 0.001 | 1.5 | 5 | 76.1 | 22.09 | 13.00 | 16.63 | 3.26 | 1.11 | 99.7% |
| 17 | 0.001 | 1.5 | 9 | 76.1 | 22.09 | 13.00 | 13.06 | 1.30 | 1.56 | 99.8% |
| 17 | 0.001 | 2 | 5 | 76.1 | 11.93 | 29.97 | 30.09 | 1.20 | 1.92 | 75.5% |
| 17 | 0.001 | 2 | 9 | 76.1 | 11.93 | 29.97 | 29.81 | 1.00 | 3.41 | 81.0% |
| 17 | 0.003 | 1.5 | 5 | 221.2 | 14.04 | 142.45 | 142.38 | 1.02 | 3.16 | 50.0% |
| 17 | 0.003 | 1.5 | 9 | 221.2 | 14.04 | 142.45 | 142.36 | 1.00 | 5.69 | 55.1% |
| 17 | 0.003 | 2 | 5 | 221.2 | 2.31 | 217.42 | 217.41 | 1.00 | 4.91 | 0.1% |
| 17 | 0.003 | 2 | 9 | 221.2 | 2.31 | 217.42 | 217.42 | 1.00 | 8.83 | 0.1% |
| 21 | 0.0005 | 1.5 | 5 | 73.4 | 28.09 | 7.63 | 15.12 | 5.01 | 1.05 | 99.9% |
| 21 | 0.0005 | 1.5 | 9 | 73.4 | 28.09 | 7.63 | 9.33 | 2.10 | 1.17 | 100.0% |
| 21 | 0.0005 | 2 | 5 | 73.4 | 20.85 | 12.96 | 16.29 | 2.91 | 1.12 | 99.9% |
| 21 | 0.0005 | 2 | 9 | 73.4 | 20.85 | 12.96 | 13.07 | 1.20 | 1.60 | 99.9% |
| 21 | 0.001 | 1.5 | 5 | 145.7 | 40.67 | 17.17 | 29.72 | 6.99 | 1.03 | 99.9% |
| 21 | 0.001 | 1.5 | 9 | 145.7 | 40.67 | 17.17 | 19.24 | 2.56 | 1.20 | 100.0% |
| 21 | 0.001 | 2 | 5 | 145.7 | 20.40 | 48.40 | 49.37 | 1.64 | 1.66 | 94.6% |
| 21 | 0.001 | 2 | 9 | 145.7 | 20.40 | 48.40 | 48.24 | 1.02 | 2.91 | 97.9% |
| 21 | 0.003 | 1.5 | 5 | 423.7 | 21.69 | 288.22 | 288.20 | 1.02 | 3.36 | 53.9% |
| 21 | 0.003 | 1.5 | 9 | 423.7 | 21.69 | 288.22 | 288.19 | 1.00 | 6.05 | 70.7% |
| 21 | 0.003 | 2 | 5 | 423.7 | 2.82 | 418.81 | 418.81 | 1.00 | 4.94 | 0.0% |
| 21 | 0.003 | 2 | 9 | 423.7 | 2.82 | 418.81 | 418.81 | 1.00 | 8.89 | 0.0% |
| 25 | 0.0005 | 1.5 | 5 | 125.2 | 47.29 | 8.97 | 25.36 | 8.99 | 1.02 | 100.0% |
| 25 | 0.0005 | 1.5 | 9 | 125.2 | 47.29 | 8.97 | 14.56 | 4.23 | 1.06 | 100.0% |
| 25 | 0.0005 | 2 | 5 | 125.2 | 34.46 | 16.33 | 25.86 | 5.62 | 1.04 | 99.9% |
| 25 | 0.0005 | 2 | 9 | 125.2 | 34.46 | 16.33 | 17.57 | 2.00 | 1.27 | 100.0% |
| 25 | 0.001 | 1.5 | 5 | 248.4 | 67.39 | 21.36 | 49.94 | 12.69 | 1.01 | 100.0% |
| 25 | 0.001 | 1.5 | 9 | 248.4 | 67.39 | 21.36 | 29.12 | 5.50 | 1.06 | 100.0% |
| 25 | 0.001 | 2 | 5 | 248.4 | 31.80 | 73.50 | 76.38 | 2.40 | 1.51 | 99.8% |
| 25 | 0.001 | 2 | 9 | 248.4 | 31.80 | 73.50 | 73.36 | 1.07 | 2.61 | 99.9% |
| 25 | 0.003 | 1.5 | 5 | 722.3 | 31.42 | 521.16 | 521.14 | 1.00 | 3.58 | 73.2% |
| 25 | 0.003 | 1.5 | 9 | 722.3 | 31.42 | 521.16 | 521.13 | 1.00 | 6.45 | 83.6% |
| 25 | 0.003 | 2 | 5 | 722.3 | 3.43 | 715.95 | 715.95 | 1.00 | 4.96 | 0.0% |
| 25 | 0.003 | 2 | 9 | 722.3 | 3.43 | 715.95 | 715.95 | 1.00 | 8.92 | 0.0% |
| 31 | 0.0005 | 1.5 | 5 | 241.3 | 89.80 | 10.79 | 48.43 | 17.70 | 1.01 | 100.0% |
| 31 | 0.0005 | 1.5 | 9 | 241.3 | 89.80 | 10.79 | 27.18 | 9.29 | 1.02 | 100.0% |
| 31 | 0.0005 | 2 | 5 | 241.3 | 64.08 | 21.33 | 48.54 | 11.96 | 1.01 | 100.0% |
| 31 | 0.0005 | 2 | 9 | 241.3 | 64.08 | 21.33 | 28.53 | 4.93 | 1.07 | 100.0% |
| 31 | 0.001 | 1.5 | 5 | 478.6 | 125.63 | 27.47 | 95.87 | 24.69 | 1.00 | 100.0% |
| 31 | 0.001 | 1.5 | 9 | 478.6 | 125.63 | 27.47 | 53.61 | 12.57 | 1.01 | 100.0% |
| 31 | 0.001 | 2 | 5 | 478.6 | 55.10 | 126.50 | 134.77 | 4.10 | 1.39 | 100.0% |
| 31 | 0.001 | 2 | 9 | 478.6 | 55.10 | 126.50 | 126.66 | 1.23 | 2.35 | 100.0% |
| 31 | 0.003 | 1.5 | 5 | 1391.6 | 50.43 | 1083.19 | 1083.19 | 1.00 | 3.88 | 85.4% |
| 31 | 0.003 | 1.5 | 9 | 1391.6 | 50.43 | 1083.19 | 1083.19 | 1.00 | 6.98 | 85.9% |
| 31 | 0.003 | 2 | 5 | 1391.6 | 4.52 | 1382.70 | 1382.70 | 1.00 | 4.97 | 0.0% |
| 31 | 0.003 | 2 | 9 | 1391.6 | 4.52 | 1382.70 | 1382.70 | 1.00 | 8.94 | 0.0% |

## Reproducing

```
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --target load_balancer_profiler -j4
build/load_balancer_profiler --d 17,21,25,31 --p 0.0005,0.001,0.003 --T 1.5,2 --k 5,9 \
    --alpha 1 --shots 100000 --warmup 1000 --seed 20260906 --out results/load_balancer
```

That writes a run directory of `shots.csv` (one row per shot, every `ticks_*` column raw, no
aggregation in the C++ at all), `agg.json` (per-cell means, maxes, percentiles, module means and
histograms, from 128-bit accumulators) and `run.log` (git hash, build flags, stim version, timer
backend and overhead, the model's assumptions, and the exact generator call and `T_int` per corpus).

Add `--verify` to cross-check the decomposed answer against the monolithic per-component path
untimed on every shot, and `--check-determinism` to diff the non-tick columns and the whole LPT
assignment across two runs.

Draw it with:

```
python benchmarks/spec_matching/plot_load_balancer.py results/load_balancer
```

Three kinds of figure come out: a per-cell histogram of `fallback` against `sparse_k` with the
`system` mean as a reference rule, a `k`-sweep small multiple per `(d, p, T)` on shared bins, and a
module breakdown per `(d, p, T)`. Add `--subtract-overhead` for the per-column timer correction
described above, and `--table` to print the summary table.

For the escalation rate `q` with a Wilson interval and a per-component `COMPLETE`/`TRUNCATED`
breakdown, use the `H`-structure profiler instead — it runs no timers and never escalates, the
trigger is the whole measurement:

```
cmake --build build --target sparse_graph_stats -j4
build/sparse_graph_stats --d 17,21,25,31 --p 0.0005,0.001,0.003 --T 1.5,2 \
    --shots 1000000 --seed 20260907 --out results/sparse_graph_stats
```

The run directory this report was written from is not in the repository — `shots.csv` alone is 4.8M
rows. It lives outside the tree, alongside it, in `lbd_apple/`.
