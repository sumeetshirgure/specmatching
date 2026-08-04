# Copyright 2026 PyReMatching contributors

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Python entry point for the two-phase decoder (§M6.3).

The C++ bindings take the detector error model as text, so that the extension does not have to
agree with a particular build of ``stim``'s own bindings. This module is the thin layer that hides
that, and the place to put horizon-unit conveniences.
"""

from typing import Any, Optional

from pyrematching._cpp_pyrematching import TwoPhaseConfig, TwoPhaseDecoder

__all__ = ["two_phase_decoder", "TwoPhaseConfig", "TwoPhaseDecoder"]


def two_phase_decoder(
    detector_error_model: Any,
    *,
    T: float = 2.0,
    T_max: Optional[float] = None,
    R: Optional[float] = None,
    config: Optional[TwoPhaseConfig] = None,
    num_distinct_weights: Optional[int] = None,
    ball_artifact_path: str = "",
) -> TwoPhaseDecoder:
    """Builds a :class:`TwoPhaseDecoder` for a ``stim.DetectorErrorModel``.

    ``T``, ``T_max`` and ``R`` are in DEM float weight units. ``T_max`` defaults to ``T`` and ``R``
    to ``2 * T_max``, which is the hard requirement of the exactness theorem (§M2.0) with no slack:
    a larger ``R`` is never wrong, a smaller one silently changes the answer.

    Pass ``config`` to set anything else; the three horizon arguments are applied on top of it.
    """
    if config is None:
        config = TwoPhaseConfig()
    config.T = T
    config.ball_T_max = T if T_max is None else T_max
    config.ball_R = 2 * config.ball_T_max if R is None else R
    if num_distinct_weights is None:
        # Whatever the extension's own default is, which is the same discretisation
        # `Matching.from_detector_error_model` uses. Naming a different one here would silently
        # give the two decoders different integer metrics.
        return TwoPhaseDecoder.from_detector_error_model(
            str(detector_error_model), config, ball_artifact_path=ball_artifact_path
        )
    return TwoPhaseDecoder.from_detector_error_model(
        str(detector_error_model),
        config,
        num_distinct_weights,
        ball_artifact_path,
    )
