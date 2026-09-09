# Latency speedups

> ## ⚠ Superseded in part — read this first
>
> **Every row in this file was produced with the lookup-table small-component resolver, and every
> row that names a `k` other than 0 is invalid: those runs decoded a *different decoder output*.**
>
> That resolver settled a connected component of `H` of size `<= k` by asking whether a perfect
> matching exists in `H` — a static, primal question. The solver decides a component by asking
> whether truncated sparse blossom *finishes by `T`* — a dynamic, dual-certified one. The two
> differ: `H` holds every edge blossom could traverse by `T` on the assumption that both endpoints
> grow for the full `T`, but blossom freezes a region the moment it is matched, so an `H` edge can
> be present and still not be reached by `T`. When that happened the solver truncated and the shot
> escalated to stock, which found the true optimum — often using an edge *outside* `H` — while the
> resolver committed `H`'s optimum, which is not the global one, and suppressed the escalation.
> Which components got which question depended on the size threshold, so the escalation rate
> depended on `k` and the output was not exact MWPM.
>
> The resolver is gone. Every component is now decided by truncated sparse blossom run on it in
> isolation, whatever its size, which is exact by the independence property: two regions interact
> only when the sum of their radii reaches the distance between their defects, every radius is
> `<= T`, so interaction needs `d_G <= 2T` — exactly the condition for an `H` edge. Defects in
> different components have no `H` edge and never interact before `T`.
>
> What this means for the tables below:
>
> * **The `k = 2` and `k = 4` rows are not measurements of any decoder that now exists.** Their
>   latencies were taken on a smaller solver input than the current path uses, and their escalation
>   rates were suppressed by the resolver's silent commits.
> * **The `escalated` column is superseded by `sparse_graph_stats`,** which measures the escalation
>   rate `q` directly, with no resolver in the path and with a per-component `COMPLETE`/`TRUNCATED`
>   breakdown behind it. Take `q` from there, not from here.
> * **The `to solver` column no longer describes anything.** Every defect goes to the solver now, so
>   that share is 1 by construction.
> * The structural columns — `components`, `largest` — describe `H` itself and are unaffected by
>   which question was asked of a component. They are the only numbers here that carry over, and
>   `sparse_graph_stats` measures them too, in more detail.
>
> **No new latency numbers were produced on this branch.** The latency profiler still charges each
> stage once per shot while the work now runs once per component, and re-reading it for that is
> separate, deferred work. This file is left in place as the record of what was measured, with the
> above stated rather than the rows quietly deleted.

Per-shot decode latency of the spec-matching decoder system against stock exact sparse blossom, over
a 48-point sweep of a rotated surface-code memory-X experiment. This file is the record of that
measurement; it replaces the per-milestone exit reports the project carried while the pipeline was
being built.

The numbers below are the **latency** read. They are not a throughput read, and the two do not agree:
on a single core running shots back to back the spec-matching path is still slower than stock exact
decode, because it does strictly more work per shot. What it buys is that most of that work is off
the critical path, and the shot can end as soon as *either* side has a usable matching.

## What was measured

Two decoders, both stock sparse blossom, timed on the same shot in a fixed order:

* **`stock`** — stock exact decode on the original detector graph `G`. No horizon, no ball tables,
  no sparsification. This is the baseline.
* **`system`** — the concurrent decoder system: `H` and `G` are solved at once and the shot ends at
  the first usable matching. On a shot the certificate rejects, that is `stock_g` — the `G` decode
  was already running, so escalation costs a lost race rather than a second decode. On every other
  shot it is `min(front end, stock_g)`.

The **front end** is `blossom + dscan + hrvst`: the solve on the sparsified graph `H`, the terminal
`max_u Y(u)` dual scan that is the certificate's own cost, and the harvest. Ball intersect, the `H`
build and the `Mwpm(H)` build are **not** charged to it — they are pipelined out. The profiler
logged their sum separately as `excluded_ns`, which its reader could add back to both series.

`speedup` throughout is `mean(stock) / mean(system)`, taken over every uncontaminated shot
**including the escalating ones**. Dropping escalating shots would price the front end at a rate no
deployment gets, so the cost of escalation sits inside the ratio rather than beside it.

### Sweep

| | |
|---|---|
| circuit | `surface_code:rotated_memory_x`, `rounds = d`, both bases decoded |
| distances `d` | 17, 21, 25, 31 |
| error rates `p` | 0.0005, 0.001, 0.003 |
| horizon `T` | 1.5, 2.0 (DEM weight units) |
| resolver cut `k` | 2, 4 |
| front end | `scan` |
| shots | 200,000 per point (48 points), 256 warmup, seed 20260907 |
| timer | `clock_thread_cputime_id`, thread-scoped |
| contamination probe | `wall_minus_thread_cpu`, 2000 ns slack |

`k` is the component-size cut: the solver saw only components of `H` larger than `k`, and sizes
`1..k` were resolved exactly off it. **That resolve is a serial pre-pass on the critical path and is
in none of these numbers** — the profiler scoped its measurement to the solver and the harvest and
left the resolve to be measured separately. The `k = 4` column is therefore an upper bound on what
the cut buys, not a delivered figure.

Shots whose wall clock ran ahead of their thread CPU clock by more than the slack were preempted
mid-decode and are dropped. Across the sweep that is 1,287,752 of 9,600,000 shots (13.4%), heavily
concentrated at `p = 0.003`, where a single shot takes long enough to be descheduled.

## Results

Speedup ranges from **1.09x** (`d=17, p=0.003, T=2`) to **6.12x** (`d=31, p=0.0005, T=1.5, k=4`).
The pattern is consistent across the grid:

* **It grows with `d`.** Stock decode on `G` scales with the whole graph; the front end scales with
  the components `H` breaks the shot into, and those do not grow the same way. At `p = 0.0005` the
  mean speedup goes 1.91 → 2.23 → 2.69 → 3.42 across `d = 17, 21, 25, 31` at `T=1.5, k=2`.
* **It shrinks with `p`.** At `p = 0.003` the sparsified graph stops being sparse: `H`'s largest
  component reaches 1082 defects at `d = 31`, over 96% of `H`'s defects are handed to the solver,
  and there is little left for the cut to remove. Speedup collapses toward 1.2–1.4.
* **`T = 1.5` beats `T = 2`, and pays for it in escalation.** The tighter horizon truncates more
  aggressively, which is where the win comes from, but escalation goes from ~0.01% to ~1% at
  `p = 0.0005`, and from 3.8% to 30.0% at `p = 0.003, d = 31`. Because escalation costs only a lost
  race and not a second decode, `T = 1.5` still wins on the mean at every point but one — `d = 31,
  p = 0.003`, where the two horizons are within 3% of each other.
* **`k = 4` beats `k = 2` nearly everywhere**, by up to 1.8x on its own (`d=31, p=0.0005, T=1.5`:
  3.42 → 6.12). The exceptions are all at `p = 0.003`, where the cut has nothing left to remove and
  the two are a wash. Subject to the caveat above that the resolve pre-pass is not charged.

### Speedup table

All latencies in microseconds; `escalated` is the share of uncontaminated shots the certificate
rejected.

| d | p | T | k | shots | stock mean | system mean | **mean x** | p50 x | p99 x | p99.9 x | escalated |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 2 | 199,758 | 6.57 | 3.43 | **1.91** | 1.96 | 1.56 | 1.54 | 0.291% |
| 17 | 0.0005 | 1.5 | 4 | 199,865 | 6.52 | 1.82 | **3.58** | 3.41 | 1.88 | 1.59 | 0.258% |
| 17 | 0.0005 | 2 | 2 | 198,814 | 7.03 | 4.96 | **1.42** | 1.42 | 1.33 | 1.49 | 0.003% |
| 17 | 0.0005 | 2 | 4 | 198,987 | 7.05 | 3.62 | **1.95** | 1.95 | 1.47 | 1.60 | 0.002% |
| 17 | 0.001 | 1.5 | 2 | 197,387 | 14.26 | 9.37 | **1.52** | 1.53 | 1.37 | 1.30 | 1.234% |
| 17 | 0.001 | 1.5 | 4 | 198,442 | 14.28 | 7.02 | **2.04** | 2.09 | 1.48 | 1.30 | 1.146% |
| 17 | 0.001 | 2 | 2 | 194,631 | 15.32 | 12.45 | **1.23** | 1.22 | 1.27 | 1.35 | 0.021% |
| 17 | 0.001 | 2 | 4 | 194,650 | 15.32 | 11.37 | **1.35** | 1.33 | 1.30 | 1.37 | 0.019% |
| 17 | 0.003 | 1.5 | 2 | 182,051 | 57.27 | 47.55 | **1.20** | 1.21 | 1.14 | 1.06 | 8.016% |
| 17 | 0.003 | 1.5 | 4 | 183,817 | 57.31 | 46.12 | **1.24** | 1.25 | 1.14 | 1.07 | 7.802% |
| 17 | 0.003 | 2 | 2 | 169,327 | 60.53 | 55.66 | **1.09** | 1.08 | 1.17 | 1.16 | 0.680% |
| 17 | 0.003 | 2 | 4 | 170,352 | 60.42 | 55.57 | **1.09** | 1.08 | 1.17 | 1.16 | 0.681% |
| 21 | 0.0005 | 1.5 | 2 | 198,437 | 13.25 | 5.95 | **2.23** | 2.27 | 1.75 | 1.49 | 0.530% |
| 21 | 0.0005 | 1.5 | 4 | 199,327 | 13.15 | 3.28 | **4.01** | 4.22 | 2.26 | 1.42 | 0.485% |
| 21 | 0.0005 | 2 | 2 | 196,247 | 15.10 | 8.77 | **1.72** | 1.73 | 1.58 | 1.64 | 0.006% |
| 21 | 0.0005 | 2 | 4 | 196,594 | 14.98 | 6.45 | **2.32** | 2.36 | 1.81 | 1.86 | 0.006% |
| 21 | 0.001 | 1.5 | 2 | 192,505 | 30.21 | 17.33 | **1.74** | 1.78 | 1.32 | 1.19 | 2.008% |
| 21 | 0.001 | 1.5 | 4 | 194,269 | 29.99 | 12.85 | **2.33** | 2.42 | 1.37 | 1.19 | 1.835% |
| 21 | 0.001 | 2 | 2 | 187,946 | 33.54 | 23.67 | **1.42** | 1.41 | 1.39 | 1.41 | 0.044% |
| 21 | 0.001 | 2 | 4 | 188,381 | 33.48 | 21.21 | **1.58** | 1.58 | 1.46 | 1.47 | 0.040% |
| 21 | 0.003 | 1.5 | 2 | 170,614 | 121.76 | 96.62 | **1.26** | 1.30 | 1.06 | 1.03 | 13.114% |
| 21 | 0.003 | 1.5 | 4 | 171,932 | 121.79 | 94.23 | **1.29** | 1.33 | 1.06 | 1.02 | 12.881% |
| 21 | 0.003 | 2 | 2 | 147,835 | 130.74 | 110.26 | **1.19** | 1.18 | 1.19 | 1.11 | 1.181% |
| 21 | 0.003 | 2 | 4 | 147,686 | 131.01 | 110.15 | **1.19** | 1.19 | 1.19 | 1.10 | 1.191% |
| 25 | 0.0005 | 1.5 | 2 | 196,683 | 25.45 | 9.45 | **2.69** | 2.77 | 1.96 | 1.27 | 0.768% |
| 25 | 0.0005 | 1.5 | 4 | 197,644 | 25.39 | 5.35 | **4.75** | 5.09 | 2.64 | 1.36 | 0.685% |
| 25 | 0.0005 | 2 | 2 | 192,742 | 28.65 | 14.30 | **2.00** | 2.01 | 1.80 | 1.78 | 0.010% |
| 25 | 0.0005 | 2 | 4 | 193,003 | 28.61 | 10.68 | **2.68** | 2.71 | 2.11 | 2.04 | 0.008% |
| 25 | 0.001 | 1.5 | 2 | 186,899 | 56.78 | 29.30 | **1.94** | 2.00 | 1.19 | 1.13 | 3.024% |
| 25 | 0.001 | 1.5 | 4 | 188,996 | 56.10 | 22.47 | **2.50** | 2.63 | 1.20 | 1.14 | 2.735% |
| 25 | 0.001 | 2 | 2 | 177,716 | 61.91 | 40.36 | **1.53** | 1.53 | 1.47 | 1.46 | 0.049% |
| 25 | 0.001 | 2 | 4 | 178,045 | 61.47 | 36.87 | **1.67** | 1.67 | 1.53 | 1.46 | 0.047% |
| 25 | 0.003 | 1.5 | 2 | 147,787 | 221.89 | 171.58 | **1.29** | 1.35 | 1.04 | 1.01 | 16.729% |
| 25 | 0.003 | 1.5 | 4 | 149,648 | 224.81 | 168.82 | **1.33** | 1.40 | 1.05 | 1.03 | 16.620% |
| 25 | 0.003 | 2 | 2 | 106,802 | 242.29 | 192.27 | **1.26** | 1.26 | 1.21 | 1.08 | 1.422% |
| 25 | 0.003 | 2 | 4 | 116,104 | 241.96 | 191.48 | **1.26** | 1.27 | 1.21 | 1.09 | 1.446% |
| 31 | 0.0005 | 1.5 | 2 | 187,955 | 60.06 | 17.54 | **3.42** | 3.57 | 1.41 | 1.19 | 1.174% |
| 31 | 0.0005 | 1.5 | 4 | 188,751 | 59.69 | 9.75 | **6.12** | 6.68 | 1.61 | 1.18 | 1.018% |
| 31 | 0.0005 | 2 | 2 | 180,463 | 66.24 | 27.04 | **2.45** | 2.46 | 2.17 | 2.05 | 0.014% |
| 31 | 0.0005 | 2 | 4 | 181,054 | 65.73 | 19.91 | **3.30** | 3.34 | 2.61 | 2.44 | 0.011% |
| 31 | 0.001 | 1.5 | 2 | 169,850 | 128.61 | 60.27 | **2.13** | 2.26 | 1.11 | 1.09 | 4.408% |
| 31 | 0.001 | 1.5 | 4 | 171,554 | 127.35 | 45.61 | **2.79** | 3.04 | 1.11 | 1.06 | 4.034% |
| 31 | 0.001 | 2 | 2 | 151,559 | 137.46 | 80.36 | **1.71** | 1.71 | 1.63 | 1.32 | 0.110% |
| 31 | 0.001 | 2 | 4 | 152,961 | 136.76 | 73.73 | **1.85** | 1.86 | 1.71 | 1.42 | 0.100% |
| 31 | 0.003 | 1.5 | 2 | 126,004 | 497.53 | 377.69 | **1.32** | 1.47 | 1.02 | 1.02 | 29.993% |
| 31 | 0.003 | 1.5 | 4 | 121,948 | 502.92 | 374.63 | **1.34** | 1.52 | 1.02 | 1.02 | 29.551% |
| 31 | 0.003 | 2 | 2 | 82,237 | 531.57 | 390.05 | **1.36** | 1.38 | 1.10 | 1.05 | 3.782% |
| 31 | 0.003 | 2 | 4 | 85,989 | 521.11 | 389.37 | **1.34** | 1.35 | 1.11 | 1.09 | 3.714% |

### Tail latencies, in microseconds

The mean is not what sizes a syndrome buffer. Blocks arrive on the hardware's clock whether or not
the previous one has been decoded, so what matters is the shape of the distribution, and the tail is
where the two series converge: at `p = 0.003` the p99 ratio is 1.02–1.21 against a mean ratio of
1.09–1.36. That is the expected shape rather than a defect — a shot in the tail is a shot with a
large component in `H`, which is exactly the shot the certificate rejects, and on a rejected shot
`system` *is* `stock_g` by construction. The front end cannot beat stock on the shots it escalates;
it can only avoid paying twice for them.

The tail ratio holds up best at low `p` and `T = 2`, where escalation is rare enough that even the
99.9th percentile is still a non-escalating shot: `d=31, p=0.0005, T=2, k=4` keeps 2.44x at p99.9.

| d | p | T | k | stock p99 | system p99 | stock p99.9 | system p99.9 |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 2 | 13.1 | 8.4 | 18.1 | 11.7 |
| 17 | 0.0005 | 1.5 | 4 | 12.9 | 6.9 | 17.4 | 11.0 |
| 17 | 0.0005 | 2 | 2 | 14.0 | 10.5 | 19.4 | 13.0 |
| 17 | 0.0005 | 2 | 4 | 14.0 | 9.5 | 19.5 | 12.2 |
| 17 | 0.001 | 1.5 | 2 | 25.5 | 18.6 | 32.7 | 25.1 |
| 17 | 0.001 | 1.5 | 4 | 25.5 | 17.2 | 32.5 | 25.0 |
| 17 | 0.001 | 2 | 2 | 27.3 | 21.5 | 35.1 | 26.0 |
| 17 | 0.001 | 2 | 4 | 27.3 | 21.0 | 34.9 | 25.4 |
| 17 | 0.003 | 1.5 | 2 | 91.3 | 80.2 | 112.0 | 106.2 |
| 17 | 0.003 | 1.5 | 4 | 91.3 | 80.1 | 112.8 | 105.6 |
| 17 | 0.003 | 2 | 2 | 97.6 | 83.6 | 120.7 | 103.9 |
| 17 | 0.003 | 2 | 4 | 97.4 | 83.6 | 119.9 | 103.4 |
| 21 | 0.0005 | 1.5 | 2 | 22.8 | 13.0 | 30.5 | 20.5 |
| 21 | 0.0005 | 1.5 | 4 | 22.4 | 9.9 | 28.2 | 19.9 |
| 21 | 0.0005 | 2 | 2 | 25.5 | 16.1 | 32.0 | 19.5 |
| 21 | 0.0005 | 2 | 4 | 25.3 | 14.0 | 32.4 | 17.4 |
| 21 | 0.001 | 1.5 | 2 | 46.0 | 34.8 | 55.1 | 46.4 |
| 21 | 0.001 | 1.5 | 4 | 45.8 | 33.4 | 54.9 | 46.0 |
| 21 | 0.001 | 2 | 2 | 50.5 | 36.3 | 60.4 | 42.8 |
| 21 | 0.001 | 2 | 4 | 50.5 | 34.6 | 60.6 | 41.3 |
| 21 | 0.003 | 1.5 | 2 | 171.2 | 161.0 | 199.1 | 194.2 |
| 21 | 0.003 | 1.5 | 4 | 171.3 | 161.2 | 198.7 | 193.8 |
| 21 | 0.003 | 2 | 2 | 184.0 | 155.1 | 212.9 | 191.5 |
| 21 | 0.003 | 2 | 4 | 184.5 | 155.3 | 213.6 | 193.8 |
| 25 | 0.0005 | 1.5 | 2 | 38.1 | 19.5 | 45.5 | 35.8 |
| 25 | 0.0005 | 1.5 | 4 | 38.2 | 14.5 | 48.2 | 35.4 |
| 25 | 0.0005 | 2 | 2 | 42.5 | 23.7 | 50.1 | 28.1 |
| 25 | 0.0005 | 2 | 4 | 42.7 | 20.2 | 50.9 | 25.0 |
| 25 | 0.001 | 1.5 | 2 | 77.9 | 65.6 | 90.4 | 79.7 |
| 25 | 0.001 | 1.5 | 4 | 77.3 | 64.4 | 89.6 | 78.9 |
| 25 | 0.001 | 2 | 2 | 84.4 | 57.3 | 97.0 | 66.4 |
| 25 | 0.001 | 2 | 4 | 83.8 | 54.9 | 96.9 | 66.3 |
| 25 | 0.003 | 1.5 | 2 | 287.2 | 276.7 | 322.6 | 318.0 |
| 25 | 0.003 | 1.5 | 4 | 297.0 | 283.0 | 341.2 | 329.8 |
| 25 | 0.003 | 2 | 2 | 314.1 | 259.7 | 350.6 | 323.9 |
| 25 | 0.003 | 2 | 4 | 312.9 | 258.5 | 351.2 | 323.7 |
| 31 | 0.0005 | 1.5 | 2 | 81.0 | 57.4 | 92.4 | 77.9 |
| 31 | 0.0005 | 1.5 | 4 | 80.1 | 49.7 | 90.9 | 76.8 |
| 31 | 0.0005 | 2 | 2 | 87.9 | 40.5 | 99.2 | 48.5 |
| 31 | 0.0005 | 2 | 4 | 87.3 | 33.4 | 97.9 | 40.2 |
| 31 | 0.001 | 1.5 | 2 | 162.1 | 146.0 | 180.4 | 165.9 |
| 31 | 0.001 | 1.5 | 4 | 159.2 | 143.3 | 173.6 | 163.7 |
| 31 | 0.001 | 2 | 2 | 170.8 | 104.8 | 184.8 | 139.8 |
| 31 | 0.001 | 2 | 4 | 170.0 | 99.2 | 185.3 | 130.6 |
| 31 | 0.003 | 1.5 | 2 | 599.0 | 589.8 | 660.3 | 645.8 |
| 31 | 0.003 | 1.5 | 4 | 614.2 | 602.2 | 672.4 | 662.3 |
| 31 | 0.003 | 2 | 2 | 640.0 | 584.2 | 699.5 | 663.4 |
| 31 | 0.003 | 2 | 4 | 634.7 | 571.1 | 709.6 | 649.4 |

### Structure of H

Same uncontaminated shots. `components` and `largest` are means over shots; `to solver` is the share
of `H`'s defects the resolver at that run's `k` had to hand over. Structural, and measured outside
every timed window — this is the explanation for the table above, not another measurement of it.

| d | p | T | k | components | largest | to solver | escalated |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 17 | 0.0005 | 1.5 | 2 | 14.91 | 6.17 | 42.3% | 0.291% |
| 17 | 0.0005 | 1.5 | 4 | 14.91 | 6.17 | 18.7% | 0.258% |
| 17 | 0.0005 | 2 | 2 | 11.31 | 9.69 | 63.9% | 0.003% |
| 17 | 0.0005 | 2 | 4 | 11.31 | 9.69 | 42.2% | 0.002% |
| 17 | 0.001 | 1.5 | 2 | 22.10 | 12.96 | 64.6% | 1.234% |
| 17 | 0.001 | 1.5 | 4 | 22.10 | 12.96 | 43.7% | 1.146% |
| 17 | 0.001 | 2 | 2 | 11.94 | 29.84 | 85.2% | 0.021% |
| 17 | 0.001 | 2 | 4 | 11.95 | 29.84 | 74.6% | 0.019% |
| 17 | 0.003 | 1.5 | 2 | 14.06 | 142.20 | 94.2% | 8.016% |
| 17 | 0.003 | 1.5 | 4 | 14.05 | 142.26 | 90.4% | 7.802% |
| 17 | 0.003 | 2 | 2 | 2.31 | 216.93 | 99.3% | 0.680% |
| 17 | 0.003 | 2 | 4 | 2.31 | 217.00 | 98.9% | 0.681% |
| 21 | 0.0005 | 1.5 | 2 | 28.08 | 7.63 | 43.2% | 0.530% |
| 21 | 0.0005 | 1.5 | 4 | 28.09 | 7.63 | 19.6% | 0.485% |
| 21 | 0.0005 | 2 | 2 | 20.83 | 12.96 | 65.5% | 0.006% |
| 21 | 0.0005 | 2 | 4 | 20.83 | 12.96 | 44.7% | 0.006% |
| 21 | 0.001 | 1.5 | 2 | 40.69 | 17.13 | 66.1% | 2.008% |
| 21 | 0.001 | 1.5 | 4 | 40.69 | 17.14 | 46.0% | 1.835% |
| 21 | 0.001 | 2 | 2 | 20.45 | 48.16 | 86.7% | 0.044% |
| 21 | 0.001 | 2 | 4 | 20.45 | 48.20 | 77.4% | 0.040% |
| 21 | 0.003 | 1.5 | 2 | 21.75 | 287.39 | 95.0% | 13.114% |
| 21 | 0.003 | 1.5 | 4 | 21.75 | 287.46 | 91.9% | 12.881% |
| 21 | 0.003 | 2 | 2 | 2.82 | 418.07 | 99.5% | 1.181% |
| 21 | 0.003 | 2 | 4 | 2.83 | 417.98 | 99.2% | 1.191% |
| 25 | 0.0005 | 1.5 | 2 | 47.25 | 8.98 | 43.9% | 0.768% |
| 25 | 0.0005 | 1.5 | 4 | 47.25 | 8.98 | 20.5% | 0.685% |
| 25 | 0.0005 | 2 | 2 | 34.43 | 16.34 | 66.7% | 0.010% |
| 25 | 0.0005 | 2 | 4 | 34.44 | 16.33 | 46.5% | 0.008% |
| 25 | 0.001 | 1.5 | 2 | 67.38 | 21.27 | 67.2% | 3.024% |
| 25 | 0.001 | 1.5 | 4 | 67.38 | 21.29 | 47.7% | 2.735% |
| 25 | 0.001 | 2 | 2 | 31.83 | 73.02 | 87.8% | 0.049% |
| 25 | 0.001 | 2 | 4 | 31.85 | 73.07 | 79.3% | 0.047% |
| 25 | 0.003 | 1.5 | 2 | 31.54 | 519.27 | 95.5% | 16.729% |
| 25 | 0.003 | 1.5 | 4 | 31.53 | 519.60 | 92.8% | 16.620% |
| 25 | 0.003 | 2 | 2 | 3.45 | 714.26 | 99.6% | 1.422% |
| 25 | 0.003 | 2 | 4 | 3.45 | 714.42 | 99.4% | 1.446% |
| 31 | 0.0005 | 1.5 | 2 | 89.76 | 10.78 | 44.8% | 1.174% |
| 31 | 0.0005 | 1.5 | 4 | 89.75 | 10.78 | 21.4% | 1.018% |
| 31 | 0.0005 | 2 | 2 | 64.08 | 21.29 | 68.1% | 0.014% |
| 31 | 0.0005 | 2 | 4 | 64.08 | 21.28 | 48.6% | 0.011% |
| 31 | 0.001 | 1.5 | 2 | 125.65 | 27.40 | 68.5% | 4.408% |
| 31 | 0.001 | 1.5 | 4 | 125.64 | 27.42 | 49.5% | 4.034% |
| 31 | 0.001 | 2 | 2 | 55.22 | 125.74 | 88.9% | 0.110% |
| 31 | 0.001 | 2 | 4 | 55.22 | 125.92 | 81.3% | 0.100% |
| 31 | 0.003 | 1.5 | 2 | 50.49 | 1081.90 | 96.1% | 29.993% |
| 31 | 0.003 | 1.5 | 4 | 50.48 | 1082.33 | 93.8% | 29.551% |
| 31 | 0.003 | 2 | 2 | 4.50 | 1386.15 | 99.7% | 3.782% |
| 31 | 0.003 | 2 | 4 | 4.50 | 1386.02 | 99.5% | 3.714% |

## Reproducing

The command below is the one that produced the rows above. **`profiler_driver` no longer exists** —
the binary and its reader `plot_latency_histograms.py` have been removed from the tree — and **`--k`
no longer exists** either, for the reason at the top of this file. The command is recorded as
provenance, not as an instruction; nothing in the current tree will run it.

```
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --target profiler_driver -j4
build/profiler_driver --distances 17,21,25,31 --error-rates 0.0005,0.001,0.003 \
    --horizons 1.5,2 --k 2,4 --modes scan --shots 200000 --seed 20260907 --out-dir results/latency
```

For the escalation rate `q` and the structure of `H`, use the `H`-structure profiler instead. It
runs no timers, keeps no decoder output and never escalates — the trigger is the whole measurement —
and it reports `q` with a Wilson interval beside a per-component `COMPLETE`/`TRUNCATED` breakdown:

```
cmake --build build --target sparse_graph_stats -j4
build/sparse_graph_stats --d 17,21,25,31 --p 0.0005,0.001,0.003 --T 1.5,2 \
    --shots 1000000 --seed 20260907 --out results/sparse_graph_stats
```

That writes `summary.csv` (one row per `(d, p, T)`), `hists.json` (the component size distribution
and the `size x status` table) and `run.log` (git hash, build flags, stim version, the exact
generator call). Add `--verify` on a smaller run to check the per-component decomposition against
the monolithic solve on every shot.

Component **size** is the only per-component statistic; the degree, hop-diameter and
weighted-diameter distributions and the boundary-structure counts were removed, from this binary and
from the library under it. Nothing about size is capped either — bin `k` counts the components of
size exactly `k`, and there is no overflow bin, so the large-component tail the escalation rate is
about is resolved rather than pooled.

Two readers take that output:

```
python benchmarks/spec_matching/plot_component_distributions.py --in results/sparse_graph_stats \
    --cells d=17,p=0.001,T=1.5
python benchmarks/spec_matching/chart_fallback_rates.py --in results/sparse_graph_stats \
    --cells d=17,p=0.001,T=1.5 d=21,p=0.001,T=1.5
```

The first draws the component size distribution of each cell on log-log axes; the second tabulates
the fraction of shots escalated, grouped by `p`, then `d`, then `T`.

The logs and figures this report was written from are not in the repository — they are ~750 MB of
per-shot CSV. They live outside the tree, alongside it, in `m1_logs_August/`.
