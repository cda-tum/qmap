"""Test cases for Clifford synthesis."""

from __future__ import annotations

import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import pytest
from qiskit import QuantumCircuit, qasm2
from qiskit.quantum_info import Clifford, PauliList

from mqt import qcec, qmap


@dataclass
class Configuration:
    """Configuration for a test case."""

    expected_minimal_gates: int
    expected_minimal_depth: int
    expected_minimal_gates_at_minimal_depth: int
    expected_minimal_two_qubit_gates: int
    expected_minimal_gates_at_minimal_two_qubit_gates: int

    description: str
    initial_tableau: str | None = None
    target_tableau: str | None = None
    initial_circuit: str | None = None
    coupling_map: str | None = None


def permute_qubits(circuit: QuantumCircuit, permutation: tuple[int, ...]) -> QuantumCircuit:
    """Return a new circuit with qubits permuted according to the given permutation."""
    permuted_circ = QuantumCircuit(circuit.num_qubits)
    qubit_map = {qubit: i for i, qubit in enumerate(circuit.qubits)}

    for d in circuit.data:
        gate = d.operation
        qubits = d.qubits
        clbits = d.clbits
        new_qubits = [circuit.qubits[permutation[qubit_map[q]]] for q in qubits]
        permuted_circ.append(gate, new_qubits, clbits)

    return permuted_circ


def convert_coupling_map(cm: str | None = None) -> list[tuple[int, int]] | None:
    """Convert a coupling map passed as a string to a CouplingMap."""
    coupling_map = None
    if cm is not None:
        pairs = cm.split(";")
        coupling_map = [(int(pair.strip("{}").split(",")[0]), int(pair.strip("{}").split(",")[1])) for pair in pairs]
    return coupling_map


def create_circuit_tests() -> list[Configuration]:
    """Create test cases for Clifford synthesis."""
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "circuits.json"
    with path.open() as f:
        circuits = json.load(f)
    return [Configuration(**c) for c in circuits]


def create_tableau_tests() -> list[Configuration]:
    """Create test cases for tableau synthesis."""
    path = Path(__file__).resolve().parent.parent / "cliffordsynthesis" / "tableaus.json"
    with path.open() as f:
        tableaus = json.load(f)
    return [Configuration(**t) for t in tableaus]


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal Clifford synthesis."""
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="gates",
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.gates == test_config.expected_minimal_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_depth(test_config: Configuration, use_maxsat: bool) -> None:
    """Test depth-optimal Clifford synthesis."""
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="depth",
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.depth == test_config.expected_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates_at_minimal_depth(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal Clifford synthesis at minimal depth."""
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="depth",
        minimize_gates_after_depth_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )
    assert results.gates == test_config.expected_minimal_gates_at_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test two-qubit gate-optimal Clifford synthesis."""
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.two_qubit_gates == test_config.expected_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_optimize_clifford_gates_at_minimal_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal Clifford synthesis at minimal two-qubit gate count."""
    circ, results = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        minimize_gates_after_two_qubit_gate_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_circuit_tests())
def test_heuristic(test_config: Configuration) -> None:
    """Test heuristic synthesis method."""
    circ, _ = qmap.optimize_clifford(
        circuit=test_config.initial_circuit,
        heuristic=True,
        split_size=10,
        target_metric="depth",
        include_destabilizers=True,
    )

    circ_opt, _ = qmap.optimize_clifford(
        circuit=test_config.initial_circuit, heuristic=False, target_metric="depth", include_destabilizers=True
    )
    assert circ.depth() >= circ_opt.depth()

    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ_opt, perm)
        if Clifford(permuted_circ) == Clifford(circ):
            equivalent = True
            break
    assert equivalent
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal tableau synthesis."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="gates",
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.gates == test_config.expected_minimal_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_depth(test_config: Configuration, use_maxsat: bool) -> None:
    """Test depth-optimal tableau synthesis."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="depth",
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.depth == test_config.expected_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates_at_minimal_depth(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal tableau synthesis at minimal depth."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="depth",
        minimize_gates_after_depth_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_depth
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test two-qubit gate-optimal tableau synthesis."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.two_qubit_gates == test_config.expected_minimal_two_qubit_gates
    print("\n", circ)


@pytest.mark.parametrize("test_config", create_tableau_tests())
@pytest.mark.parametrize("use_maxsat", [True, False], ids=["maxsat", "binary_search"])
def test_synthesize_clifford_gates_at_minimal_two_qubit_gates(test_config: Configuration, use_maxsat: bool) -> None:
    """Test gate-optimal tableau synthesis at minimal two-qubit gate count."""
    circ, results = qmap.synthesize_clifford(
        target_tableau=test_config.target_tableau,
        initial_tableau=test_config.initial_tableau,
        use_maxsat=use_maxsat,
        target_metric="two_qubit_gates",
        try_higher_gate_limit_for_two_qubit_gate_optimization=True,
        minimize_gates_after_two_qubit_gate_optimization=True,
        coupling_map=convert_coupling_map(test_config.coupling_map),
    )

    assert results.gates == test_config.expected_minimal_gates_at_minimal_two_qubit_gates
    print("\n", circ)


# The following tests merely check that all kinds of conversions work as expected.
# They only check for the correctness of the result for a simple circuit.


@pytest.fixture
def bell_circuit() -> QuantumCircuit:
    """Return a Bell state circuit."""
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.cx(0, 1)
    return circ


def test_optimize_quantum_computation(bell_circuit: QuantumCircuit) -> None:
    """Test that we can optimize an MQT QuantumComputation."""
    qc = qmap.QuantumComputation.from_qiskit(bell_circuit)
    circ, _ = qmap.optimize_clifford(circuit=qc)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_optimize_from_qasm_file(bell_circuit: QuantumCircuit) -> None:
    """Test that we can optimize from a QASM file."""
    qasm2.dump(bell_circuit, Path("bell.qasm"))
    circ, _ = qmap.optimize_clifford(circuit="bell.qasm")
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_optimize_qiskit_circuit(bell_circuit: QuantumCircuit) -> None:
    """Test that we can optimize a Qiskit QuantumCircuit."""
    circ, _ = qmap.optimize_clifford(circuit=bell_circuit)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_optimize_with_initial_tableau(bell_circuit: QuantumCircuit) -> None:
    """Test that we can optimize a circuit with an initial tableau."""
    circ, _ = qmap.optimize_clifford(circuit=bell_circuit, initial_tableau=qmap.Tableau(bell_circuit.num_qubits))
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_optimize_with_initial_tableau_with_mapping(bell_circuit: QuantumCircuit) -> None:
    """Test that we can optimize a circuit with an initial tableau."""
    coupling_map = [(0, 1), (1, 0)]
    circ, _ = qmap.optimize_clifford(
        circuit=bell_circuit, initial_tableau=qmap.Tableau(bell_circuit.num_qubits), coupling_map=coupling_map
    )
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_from_tableau(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from an MQT Tableau."""
    tableau = qmap.Tableau("['XX', 'ZZ']")
    circ, _ = qmap.synthesize_clifford(target_tableau=tableau)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_from_qiskit_clifford(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from a Qiskit Clifford."""
    cliff = Clifford(bell_circuit)
    circ, _ = qmap.synthesize_clifford(target_tableau=cliff)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_from_qiskit_pauli_list(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from a Qiskit PauliList."""
    pauli_list = PauliList(["XX", "ZZ"])
    circ, _ = qmap.synthesize_clifford(target_tableau=pauli_list)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_from_string(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from a String."""
    pauli_str = "[XX,ZZ]"
    circ, _ = qmap.synthesize_clifford(target_tableau=pauli_str)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_with_coupling_from_tableau(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from an MQT Tableau."""
    tableau = qmap.Tableau("['XX', 'ZZ']")
    coupling_map = [(0, 1), (1, 0)]
    circ, _ = qmap.synthesize_clifford(target_tableau=tableau, coupling_map=coupling_map)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_with_coupling_with_initial_tableau_from_tableau(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from an MQT Tableau."""
    tableau = qmap.Tableau("['XX', 'ZZ']")
    init_tableau = qmap.Tableau("['ZI','IZ']")
    coupling_map = [(0, 1), (1, 0)]
    circ, _ = qmap.synthesize_clifford(target_tableau=tableau, initial_tableau=init_tableau, coupling_map=coupling_map)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_synthesize_with_coupling_with_initial_tableau_from_qc(bell_circuit: QuantumCircuit) -> None:
    """Test that we can synthesize a circuit from an MQT Tableau."""
    init_tableau = qmap.Tableau("['ZI','IZ']")
    coupling_map = [(0, 1), (1, 0)]
    circ, _ = qmap.optimize_clifford(circuit=bell_circuit, initial_tableau=init_tableau, coupling_map=coupling_map)
    num_qubits = circ.num_qubits
    qubit_permutations = list(itertools.permutations(range(num_qubits)))
    equivalent = False
    for perm in qubit_permutations:
        permuted_circ = permute_qubits(circ, perm)
        if qcec.verify(permuted_circ, bell_circuit).considered_equivalent():
            equivalent = True
            break
    assert equivalent


def test_invalid_kwarg_to_synthesis() -> None:
    """Test that we raise an error if we pass an invalid kwarg to synthesis."""
    with pytest.raises(ValueError, match="Invalid keyword argument"):
        qmap.synthesize_clifford(target_tableau=qmap.Tableau("Z"), invalid_kwarg=True)
