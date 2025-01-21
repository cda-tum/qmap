"""Test the state preparation for zoned neutral atom architectures."""

from __future__ import annotations

from pathlib import Path

import pytest

from mqt.core import load
from mqt.qmap.na import NAStatePreparationSolver, generate_code, get_ops_for_solver

# get root directory of the project
circ_dir = Path(__file__).resolve().parent.parent.parent / "na/nasp/circuits"


@pytest.fixture
def solver() -> NAStatePreparationSolver:
    """Return a NAStatePreparationSolver instance with both sided storage zone."""
    return NAStatePreparationSolver(3, 7, 2, 3, 2, 2, 2, 2, 2, 4)


@pytest.mark.parametrize(
    ("circuit_filename", "n_qubits"),
    [
        ("steane.qasm", 7),
        ("surface_3.qasm", 9),
    ],
)
def test_na_state_prep_sat(solver: NAStatePreparationSolver, circuit_filename: str, n_qubits: int) -> None:
    """Test the state preparation for the zoned neutral atom architecture."""
    qc = load(circ_dir / circuit_filename)
    ops = get_ops_for_solver(qc, "z", 1)
    assert ops is not None
    assert len(ops) > 0
    result = solver.solve(ops, n_qubits, 4, None, False, True)
    assert result is not None
    assert result.yaml().startswith("sat: true")
    code = generate_code(qc, result, 1, 10, 24)
    assert code is not None
    assert len(code.splitlines()) >= 7  # for each stage there is at least one line


@pytest.mark.parametrize(
    ("circuit_filename", "n_qubits"),
    [
        ("shor.qasm", 9),
    ],
)
def test_na_state_prep_unsat(solver: NAStatePreparationSolver, circuit_filename: str, n_qubits: int) -> None:
    """Test the state preparation for the zoned neutral atom architecture."""
    qc = load(circ_dir / circuit_filename)
    ops = get_ops_for_solver(qc, "z", 1)
    assert ops is not None
    assert len(ops) > 0
    result = solver.solve(ops, n_qubits, 4, None, False, True)
    assert result is not None
    assert result.yaml().startswith("sat: false")
