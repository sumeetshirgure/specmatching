# SpecMatching

SpecMatching is a Python/C++ library for accelerating the PyMatching decoder with speculative pre-decoding.

![SpecBlossom](data/specblossom.png)

## The spec-matching decoder

Dual - truncated sparse blossom, executed on a precomputed ball graph over the shot's defects, with stock
exact decode as the residual fallback. **The output is exact MWPM on every shot** — the fallback is
not an approximation, it is the inherited decoder run in full on the shots the truncated phase did
not finish.

`T` is the truncation horizon; the ball tables are sized from it and take `R >= 2 * T` of radius,
which is the memory cost. Lowering `T` raises the fraction of shots that fall back — measured at
`5e-6` to `1.3e-4` at `T = 2` and one to two orders higher at `T = 1.5`.

`specmatching.spec_matching_decoder(dem, T=2.0, stock_on_h=True)` selects an alternative front end:
**stock** (untruncated) sparse blossom on the ball graph, kept when its terminal duals certify the
shot globally optimal and escalated otherwise. Same exact output, same escalation rate (measured
identical over 10⁶ shots), about 3% slower, and a cleaner correctness argument — it is checked by
direct equality against stock exact decode rather than against a truncated intermediate state.

**Latency, not throughput.** Run as a system — the sparsified graph `H` and the original graph `G`
solved concurrently, the shot ending at the first usable matching — per-shot latency is **1.1x to
6.1x** faster than stock exact decode across `d = 17,21,25,31` and `p = 0.0005,0.001,0.003`.
On a single core running shots back to back it is still *slower* than stock, because it does strictly
more work per shot; what it buys is that most of that work is off the critical path.
The full 48-point sweep, the tail percentiles and the method are in [`docs/latency_speedups.md`](docs/latency_speedups.md).

Currently this is not distributed as a software tool because ball graph pre-computation
is memory bound on CPU.

## Attribution

The paper for the speculation-accelerated version is coming out soon.

When using SpecMatching please also cite the original [paper](https://arxiv.org/abs/2303.15933) on the sparse blossom
as `specmatching` is build on `pymatching`:

```
@article{Higgott2025sparseblossom,
  doi = {10.22331/q-2025-01-20-1600},
  url = {https://doi.org/10.22331/q-2025-01-20-1600},
  title = {Sparse {B}lossom: correcting a million errors per core second with minimum-weight matching},
  author = {Higgott, Oscar and Gidney, Craig},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {9},
  pages = {1600},
  month = jan,
  year = {2025}
}
```
