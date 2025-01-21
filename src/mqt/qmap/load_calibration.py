"""Module for loading calibration data."""

from __future__ import annotations

from typing import TYPE_CHECKING

from qiskit.providers.models import BackendProperties
from qiskit.transpiler.target import Target

if TYPE_CHECKING:
    from .pyqmap import Architecture


def load_calibration(architecture: Architecture, calibration: str | Target | BackendProperties | None = None) -> None:
    """Load a calibration from a string, BackendProperties, or Target.

    Args:
        architecture: The architecture to load the calibration into.
        calibration: The calibration to load.
    """
    if calibration is None:
        return

    if isinstance(calibration, str):
        architecture.load_properties(calibration)
    elif isinstance(calibration, BackendProperties):
        from mqt.qmap.plugins.qiskit import import_backend_properties

        architecture.load_properties(import_backend_properties(calibration))
    elif isinstance(calibration, Target):
        from mqt.qmap.plugins.qiskit import import_target

        architecture.load_properties(import_target(calibration))
    else:  # pragma: no cover
        msg = f"Calibration type {type(calibration)} not supported."
        raise TypeError(msg)
