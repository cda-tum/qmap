"""MQT QMAP library.

This file is part of the MQT QMAP library released under the MIT license.
See README.md or go to https://github.com/cda-tum/qmap for more information.
"""

from __future__ import annotations

import sys

# under Windows, make sure to add the appropriate DLL directory to the PATH
if sys.platform == "win32":

    def _dll_patch() -> None:
        """Add the DLL directory to the PATH."""
        import os
        import sysconfig
        from pathlib import Path

        bin_dir = Path(sysconfig.get_paths()["purelib"]) / "mqt" / "core" / "bin"
        os.add_dll_directory(str(bin_dir))

        if sys.version_info >= (3, 9, 0) and "Z3_ROOT" in os.environ:
            lib_path = Path(os.environ["Z3_ROOT"]) / "lib"
            if lib_path.exists():
                os.add_dll_directory(str(lib_path))
            bin_path = Path(os.environ["Z3_ROOT"]) / "bin"
            if bin_path.exists():
                os.add_dll_directory(str(bin_path))

    _dll_patch()
    del _dll_patch

from ._version import version as __version__
from .clifford_synthesis import optimize_clifford, synthesize_clifford
from .compile import compile  # noqa: A004
from .pyqmap import (
    Arch,
    Architecture,
    CliffordSynthesizer,
    CommanderGrouping,
    Configuration,
    EarlyTermination,
    Encoding,
    Heuristic,
    HybridMapperParameters,
    HybridNAMapper,
    InitialCircuitMapping,
    InitialCoordinateMapping,
    InitialLayout,
    Layering,
    LookaheadHeuristic,
    MappingResults,
    Method,
    NeutralAtomHybridArchitecture,
    SwapReduction,
    SynthesisConfiguration,
    SynthesisResults,
    Tableau,
    TargetMetric,
    Verbosity,
)
from .subarchitectures import SubarchitectureOrder

__all__ = [
    "Arch",
    "Architecture",
    "CliffordSynthesizer",
    "CommanderGrouping",
    "Configuration",
    "EarlyTermination",
    "Encoding",
    "Heuristic",
    "HybridMapperParameters",
    "HybridNAMapper",
    "InitialCircuitMapping",
    "InitialCoordinateMapping",
    "InitialLayout",
    "Layering",
    "LookaheadHeuristic",
    "MappingResults",
    "Method",
    "NeutralAtomHybridArchitecture",
    "SubarchitectureOrder",
    "SwapReduction",
    "SynthesisConfiguration",
    "SynthesisResults",
    "Tableau",
    "TargetMetric",
    "Verbosity",
    "__version__",
    "compile",
    "optimize_clifford",
    "synthesize_clifford",
]
