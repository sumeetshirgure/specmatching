# Copyright 2026 SpecMatching contributors

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""§M6 exit checkpoint: the bindings are usable from python end to end, artifact save/load
included, and `summarize()` reconciles the amortised mean against the measured mean."""

import os
import tempfile

import numpy as np
import pytest
import stim

import specmatching


def _circuit(distance=5, rounds=5, noise=0.01):
    return stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=noise,
        after_reset_flip_probability=noise,
        before_measure_flip_probability=noise,
    )


def _sample(circuit, shots, seed=20260812):
    sampler = circuit.compile_detector_sampler(seed=seed)
    dets, obs = sampler.sample(shots, separate_observables=True)
    return [np.flatnonzero(row).astype(np.uint64) for row in dets], obs


@pytest.fixture(scope="module")
def corpus():
    circuit = _circuit()
    dem = circuit.detector_error_model(decompose_errors=True)
    shots, _ = _sample(circuit, 200)
    return dem, shots


def _stock_predictions(dem, shots):
    matching = specmatching.Matching.from_detector_error_model(dem)
    out = []
    for shot in shots:
        obs, weight = matching.decode(
            np.isin(np.arange(dem.num_detectors), shot).astype(np.uint8), return_weight=True
        )
        out.append((np.asarray(obs, dtype=np.uint8), weight))
    return out


def test_decode_is_exact_against_stock(corpus):
    dem, shots = corpus
    decoder = specmatching.spec_matching_decoder(dem, T=2.0)
    expected = _stock_predictions(dem, shots)
    for shot, (want_obs, want_weight) in zip(shots, expected):
        got_obs, got_weight = decoder.decode_to_obs(shot)
        np.testing.assert_array_equal(got_obs, want_obs)
        # The two weights are the same integer in the discretised metric, reported as a float.
        assert got_weight == pytest.approx(want_weight, rel=1e-9, abs=1e-9)


def test_decode_batch_and_profile_are_row_aligned(corpus):
    dem, shots = corpus
    decoder = specmatching.spec_matching_decoder(dem, T=2.0)
    obs, weights, profile = decoder.decode_batch(shots, profile=True)
    assert obs.shape == (len(shots), decoder.num_observables)
    assert weights.shape == (len(shots),)
    assert profile is not None
    for column in profile.values():
        assert len(column) == len(shots)

    for i, shot in enumerate(shots):
        single_obs, single_weight = decoder.decode_to_obs(shot)
        np.testing.assert_array_equal(obs[i], single_obs)
        assert weights[i] == single_weight

    # Invariant 12, visible from python: escalation fires exactly when the timeline truncated.
    np.testing.assert_array_equal(profile["escalated"], profile["truncated"])


def test_summarize_reconciles_the_amortised_mean(corpus):
    dem, shots = corpus
    config = specmatching.SpecMatchingConfig()
    config.measure_exact_reference = True
    decoder = specmatching.spec_matching_decoder(dem, T=2.0, config=config)
    decoder.decode_batch(shots, profile=True)

    summary = specmatching.summarize(decoder.get_aggregate_stats())
    assert summary["shots"] == len(shots)
    # `q` is meaningless without its denominator: a zero below 1/shots is a resolution floor.
    assert summary["q"] == summary["shots_escalated"] / summary["shots"]
    assert summary["measured_mean_ns"] > 0
    # The M6 exit checkpoint's reconciliation. The gap is unattributed cost — the parts the profile
    # does not name — and on this path it is small but not zero.
    assert summary["amortisation_gap"] < 0.5
    assert summary["speedup_vs_stock"] > 0
    assert summary["p9999_total_ns"] >= summary["p999_total_ns"] >= summary["p99_total_ns"]


def test_edges_flavor_is_a_valid_correction(corpus):
    dem, shots = corpus
    config = specmatching.SpecMatchingConfig()
    config.edges_flavor = True
    config.store_paths = True
    decoder = specmatching.spec_matching_decoder(dem, T=2.0, config=config)

    for shot in shots[:50]:
        edges, _ = decoder.decode_to_edges(shot)
        parity = np.zeros(dem.num_detectors, dtype=np.uint8)
        for u, v in edges:
            parity[u] ^= 1
            if v >= 0:
                parity[v] ^= 1
        np.testing.assert_array_equal(np.flatnonzero(parity), np.sort(shot))


def test_ball_artifact_round_trips(corpus):
    dem, shots = corpus
    decoder = specmatching.spec_matching_decoder(dem, T=2.0)
    stats = decoder.ball_stats()
    assert stats["bytes_total"] > 0
    assert stats["num_nodes"] == dem.num_detectors

    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "balls.artifact")
        decoder.save_ball_artifact(path)
        loaded = specmatching.spec_matching_decoder(dem, T=2.0, ball_artifact_path=path)
        for shot in shots[:50]:
            np.testing.assert_array_equal(loaded.decode_to_obs(shot)[0], decoder.decode_to_obs(shot)[0])


def test_stock_on_h_is_exact_against_stock(corpus):
    """§M7 from python: the certificate front end is exact MWPM on every shot, and reports itself."""
    dem, shots = corpus
    decoder = specmatching.spec_matching_decoder(dem, T=2.0, stock_on_h=True)
    expected = _stock_predictions(dem, shots)
    for shot, (want_obs, want_weight) in zip(shots, expected):
        got_obs, got_weight = decoder.decode_to_obs(shot)
        np.testing.assert_array_equal(got_obs, want_obs)
        assert got_weight == pytest.approx(want_weight, rel=1e-9, abs=1e-9)

    _, _, profile = decoder.decode_batch(shots, profile=True)
    # Exactly one of the two branches, on every shot.
    np.testing.assert_array_equal(profile["certified"] == 1, ~profile["escalated"])
    # A shot that escalated because `H` had no perfect matching has no dual to report.
    assert np.all(profile["max_dual_at_completion"][profile["h_no_perfect_matching"] == 1] == 0)


def test_stock_on_h_reports_q_this_against_q_current(corpus):
    """§M7.7: the landed scheme's decision replayed on identical shots, so the reduction is measured.

    `q_current_on_same_corpus` is `None` rather than `0` when the replay was not run — "not measured"
    and "measured, zero" are different claims, and at these rates the second is a resolution floor.
    """
    dem, shots = corpus
    config = specmatching.SpecMatchingConfig()
    config.stock_on_h = True
    config.measure_truncated_reference = True
    decoder = specmatching.spec_matching_decoder(dem, T=1.0, config=config)
    decoder.decode_batch(shots, profile=True)

    summary = specmatching.summarize(decoder.get_aggregate_stats())
    assert summary["q_this"] == summary["q"]
    assert summary["shots_with_truncated_reference"] == len(shots)
    assert summary["q_current_on_same_corpus"] is not None
    # §M7.0's corollary, on identical shots.
    assert summary["q"] <= summary["q_current_on_same_corpus"]

    without_replay = specmatching.spec_matching_decoder(dem, T=1.0, stock_on_h=True)
    without_replay.decode_batch(shots, profile=True)
    assert specmatching.summarize(without_replay.get_aggregate_stats())["q_current_on_same_corpus"] is None


def test_stock_on_h_rejects_the_m1_harvest_oracle(corpus):
    """The M1 oracle has no truncated intermediate state to reproduce here; M7's oracle is stock-on-G."""
    dem, _ = corpus
    config = specmatching.SpecMatchingConfig()
    config.stock_on_h = True
    config.verify_against_g = True
    with pytest.raises(ValueError):
        specmatching.spec_matching_decoder(dem, T=2.0, config=config)


def test_unbounded_horizon_requires_the_oracle_front_end(corpus):
    dem, _ = corpus
    config = specmatching.SpecMatchingConfig()
    config.unbounded_horizon = True
    with pytest.raises(ValueError):
        specmatching.spec_matching_decoder(dem, T=2.0, config=config)
