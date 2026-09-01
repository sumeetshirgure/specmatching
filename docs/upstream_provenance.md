# Upstream provenance

`src/specmatching/`, excluding `spec_matching/`, is [PyMatching](https://github.com/oscarhiggott/PyMatching)
vendored in and renamed wholesale, per §M0 of the implementation plan. This file is the difference
between a fork and a laundering: it says what was taken, from where, and how.

| | |
|---|---|
| Upstream | `https://github.com/oscarhiggott/PyMatching` |
| Upstream version | **v2.4.0** — the version `src/specmatching/_version.py` carried before the reset below, and the version §M0.0's ground-truth table was verified against |
| Upstream commit SHA | **not recorded.** See "the gap" below |
| Vendoring commit | `66d374669267f01fd53f7993e9ac45535e1c7499`, *"Fork" pymatching and rename to specmatching*, tagged `vendor-base` |
| Vendoring command | `python tools/vendor_pymatching.py`, verified by `tools/check_vendor.sh` |
| `extern/pybind11` | pinned at commit `59565095fafae453b15d2a9e8a66d8ca4758b6d5` (`v3.0.2-79-g59565095`) |
| Licence | Apache-2.0. `LICENSE` is upstream's, verbatim, with its copyright line intact. `NOTICE` names PyMatching, its copyright holder and the upstream URL, and identifies `src/specmatching/sparse_blossom/` as derived work |
| Package version | reset to `0.0.1` (§M0.5). Inheriting upstream's `2.4.0` on PyPI would be confusing and would look like a hijack |

Every subsequent edit to a vendored file has a row in [`upstream_edits.md`](upstream_edits.md).
That file, not this one, is the authoritative list of divergence.

## The gap, stated rather than papered over

§M0.1 step 2 asks for the upstream commit SHA to be recorded at clone time, on the grounds that it
is "the only way anyone will ever diff you against upstream again". **It was not recorded**, and it
cannot be reconstructed from this repository: the vendoring commit's message does not carry it, no
upstream remote exists, and the vendoring script does not stamp it.

What survives is the *version*, v2.4.0, which pins the diff to a release tag rather than to a
commit. That is enough to re-derive the vendored tree — `tools/vendor_pymatching.py` is idempotent
by design, so cloning upstream at `v2.4.0` and re-running it should reproduce
`66d3746`'s tree — and it is not enough to attribute any post-tag upstream commit that may have been
included. Anyone doing that diff should start by checking whether the re-derived tree matches, and
record the answer here.

This is written down rather than quietly fixed with a plausible SHA, because a provenance file whose
SHA is a guess is worse than one that says it does not know.

## `check_vendor.sh` gate status

| Gate | Status |
|---|---|
| 1. no stray `pymatching` identifiers | Passes over the tracked tree. The script also scans `design/`, which is untracked and gitignored and is the *design document* — it discusses upstream by name throughout, and every hit there is intentional prose |
| 2. upstream URLs survived the rename | Skipped unless `--upstream <path-to-a-pymatching-clone>` is passed. Run it with a clone at `v2.4.0` before any release |
| 3. builds | Passes |
| 4. imports and decodes | Passes |
