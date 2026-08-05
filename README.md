# PyReMatching

PyReMatching is a Python/C++ library for accelerating the PyMatching decoder.

## The two-phase decoder

Truncated sparse blossom, executed on a precomputed ball graph over the shot's defects, with stock
exact decode as the residual fallback. **The output is exact MWPM on every shot** — the fallback is
not an approximation, it is the inherited decoder run in full on the shots the truncated phase did
not finish.

```python
import pyrematching
import stim

circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=11, rounds=11,
                                 after_clifford_depolarization=1e-3)
dem = circuit.detector_error_model(decompose_errors=True)

decoder = pyrematching.two_phase_decoder(dem, T=2.0)   # T is in DEM weight units
observables, weight = decoder.decode_to_obs(detection_event_indices)

observables, weights, profile = decoder.decode_batch(shots, profile=True)
print(pyrematching.summarize(decoder.get_aggregate_stats()))
```

`T` is the truncation horizon; the ball tables are sized from it and take `R >= 2 * T` of radius,
which is the memory cost. Lowering `T` raises the fraction of shots that fall back — measured at
`5e-6` to `1.3e-4` at `T = 2` and one to two orders higher at `T = 1.5`.

`pyrematching.two_phase_decoder(dem, T=2.0, stock_on_h=True)` selects an alternative front end:
**stock** (untruncated) sparse blossom on the ball graph, kept when its terminal duals certify the
shot globally optimal and escalated otherwise. Same exact output, same escalation rate (measured
identical over 10⁶ shots), about 3% slower, and a cleaner correctness argument — it is checked by
direct equality against stock exact decode rather than against a truncated intermediate state. See
`docs/two_phase_m7_exit.md`.

**It is not currently faster than the inherited decoder on a CPU**: measured at 0.28–0.42x stock
exact decode throughput. See `docs/two_phase_m6_exit.md` for why, and for the critical-path and
local-memory measurements the design is actually argued on. The exit reports for each milestone are
in `docs/two_phase_m*_exit.md`.


## Attribution

When using PyReMatching please cite the original [paper](https://arxiv.org/abs/2303.15933) on the sparse blossom algorithm (implemented in version 2):

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


The paper for the accelerated version is coming soon.
