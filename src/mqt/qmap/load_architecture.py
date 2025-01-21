"""Module for loading architectures."""

from __future__ import annotations

from qiskit.providers import Backend

from .pyqmap import Arch, Architecture

__all__ = [
    "load_architecture",
]


def __dir__() -> list[str]:
    return __all__


def load_architecture(arch: str | Arch | Architecture | Backend | None = None) -> Architecture:
    """Load an architecture from a string, Arch, Architecture, or Backend. If None is passed, no architecture is loaded.

    Args:
        arch: The architecture to load.

    Returns:
        The loaded architecture.
    """
    architecture = Architecture()

    if arch is not None:
        if isinstance(arch, str):
            try:
                architecture.load_coupling_map(Arch(arch))
            except ValueError:
                architecture.load_coupling_map(arch)
        elif isinstance(arch, Arch):
            architecture.load_coupling_map(arch)
        elif isinstance(arch, Architecture):
            architecture = arch
        elif isinstance(arch, Backend):
            from mqt.qmap.plugins.qiskit import import_backend

            architecture = import_backend(arch)
        else:  # pragma: no cover
            msg = f"Architecture type {type(arch)} not supported."
            raise TypeError(msg)

    return architecture
