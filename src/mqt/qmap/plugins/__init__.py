"""Plugins for the QMAP package."""

from __future__ import annotations

from . import qiskit

__all__ = [
    "qiskit",
]


def __dir__() -> list[str]:
    return __all__
