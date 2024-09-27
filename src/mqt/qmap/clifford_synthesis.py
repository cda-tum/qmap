"""Main entry point for the Clifford synthesis module."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from qiskit import QuantumCircuit, qasm3
from qiskit.quantum_info import Clifford, PauliList
from qiskit.transpiler.layout import TranspileLayout

from .compile import extract_initial_layout_from_qasm
from .pyqmap import (
    CliffordSynthesizer,
    QuantumComputation,
    SynthesisConfiguration,
    SynthesisResults,
    Tableau,
)

if TYPE_CHECKING:
    from collections.abc import Iterable


def _import_circuit(circuit: str | QuantumCircuit | QuantumComputation) -> QuantumComputation:
    """Import a circuit from a string, a QuantumCircuit, or a QuantumComputation."""
    if isinstance(circuit, QuantumCircuit):
        return QuantumComputation.from_qiskit(circuit)
    if isinstance(circuit, str):
        if circuit.endswith(".qasm"):
            return QuantumComputation.from_file(circuit)
        return QuantumComputation.from_qasm_str(circuit)
    return circuit


def _reverse_paulis(paulis: list[str]) -> list[str]:
    return [s[0] + s[:0:-1] if s[0] in "+-" else s[::-1] for s in paulis]


def _import_tableau(tableau: str | Clifford | PauliList | Tableau, include_destabilizers: bool = False) -> Tableau:
    """Import a tableau from a string, a Clifford, a PauliList, or a Tableau."""
    if isinstance(tableau, Clifford):
        mode = "B" if include_destabilizers else "S"
        return Tableau(str(_reverse_paulis(tableau.to_labels(mode=mode))))
    if isinstance(tableau, PauliList):
        return Tableau(str(_reverse_paulis(tableau.to_labels())))
    if isinstance(tableau, str):
        return Tableau(tableau)
    return tableau


def _config_from_kwargs(kwargs: dict[str, Any]) -> SynthesisConfiguration:
    """Create a :class:`SynthesisConfiguration` from keyword arguments."""
    config = SynthesisConfiguration()
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)
        else:
            msg = f"Invalid keyword argument: {key}"
            raise ValueError(msg)

    if not config.solver_parameters:
        config.solver_parameters = {}
        if config.use_maxsat:
            config.solver_parameters["pb.compile_equality"] = True
            config.solver_parameters["maxres.hill_climb"] = True
            config.solver_parameters["maxres.pivot_on_correction_set"] = False
        else:
            config.solver_parameters["bca"] = True
            config.solver_parameters["restart.emafastglue"] = 0.05
            config.solver_parameters["restart.emaslowglue"] = 1e-6
            config.solver_parameters["restart.margin"] = 1.07
            config.solver_parameters["rephase.base"] = 3000
            config.solver_parameters["search.sat.conflicts"] = 100

    return config


def _circuit_from_qasm(qasm: str) -> QuantumCircuit:
    """Create a proper :class:`qiskit.QuantumCircuit` from a QASM string (including layout information)."""
    circ = qasm3.loads(qasm)
    layout = extract_initial_layout_from_qasm(qasm, circ.qregs)

    circ._layout = TranspileLayout(  # noqa: SLF001
        initial_layout=layout, input_qubit_mapping=layout.get_virtual_bits()
    )

    return circ


def synthesize_clifford(
    target_tableau: str | Clifford | PauliList | Tableau,
    initial_tableau: str | Clifford | PauliList | Tableau | None = None,
    include_destabilizers: bool = False,
    coupling_map: Iterable[tuple[int, int]] | None = None,
    **kwargs: Any,  # noqa: ANN401
) -> tuple[QuantumCircuit, SynthesisResults]:
    """Synthesize a Clifford circuit from a given tableau starting from an (optional) initial tableau.

    Args:
        target_tableau:
            The target tableau to synthesize.
            If a string is given, it is interpreted as a semicolon separated binary matrix or a list of Pauli strings. The Pauli strings follow the same format as in `Stim <https://github.com/quantumlib/Stim>`_.
            If a :class:`Clifford` or a :class:`PauliList` is given, it is converted to a :class:`Tableau`.
            If a :class:`Tableau` is given, it is used directly.
        initial_tableau:
            The initial tableau to start from.
            If a string is given, it is interpreted as a semicolon separated binary matrix or a list of Pauli strings.
            If a :class:`Clifford` or a :class:`PauliList` is given, it is converted to a :class:`Tableau`.
            If a :class:`Tableau` is given, it is used directly.
            If no initial tableau is given, the synthesis starts from the identity tableau.
        include_destabilizers:
            Flag to set whether destabilizers should be considered in the synthesis
        coupling_map:
             Iterable of tuples of connected qubits of the architecture to synthesize for.
        **kwargs:
            Additional keyword arguments to configure the synthesis.
            See :class:`SynthesisConfiguration` for a list of available options.

    Returns:
        A tuple containing the synthesized circuit and the synthesis results.
    """
    config = _config_from_kwargs(kwargs)

    tableau = _import_tableau(target_tableau, include_destabilizers)
    if initial_tableau is not None:
        if coupling_map is not None:
            synthesizer = CliffordSynthesizer(_import_tableau(initial_tableau), tableau, set(coupling_map))
        else:
            synthesizer = CliffordSynthesizer(_import_tableau(initial_tableau), tableau)
    elif coupling_map is not None:
        synthesizer = CliffordSynthesizer(tableau, set(coupling_map))
    else:
        synthesizer = CliffordSynthesizer(tableau)

    synthesizer.synthesize(config)

    results = synthesizer.results
    circ = _circuit_from_qasm(results.circuit)

    return circ, results


def optimize_clifford(
    circuit: str | QuantumCircuit | QuantumComputation,
    initial_tableau: str | Clifford | PauliList | Tableau | None = None,
    include_destabilizers: bool = False,
    coupling_map: Iterable[tuple[int, int]] | None = None,
    **kwargs: Any,  # noqa: ANN401
) -> tuple[QuantumCircuit, SynthesisResults]:
    """Optimize a Clifford circuit starting from an (optional) initial tableau.

    Args:
        circuit:
            The circuit to optimize.
            If a string is given, it is interpreted as a QASM string or a filename.
            If a :class:`QuantumCircuit` is given, it is converted to a :class:`QuantumComputation`.
            If a :class:`QuantumComputation` is given, it is used as is.
        initial_tableau:
            The initial tableau to start from.
            If a string is given, it is interpreted as a semicolon separated binary matrix or a list of Pauli strings.
            If a :class:`Clifford` is given or a :class:`PauliList` is given, it is converted to a Tableau.
            If a :class:`Tableau` is given, it is used directly.
            If no initial tableau is given, the synthesis starts from the identity tableau.
        include_destabilizers:
            Flag to set whether destabilizers should be considered in the synthesis
        coupling_map:
            Iterable of tuples of connected qubits of the architecture to synthesize for.
        **kwargs:
            Additional keyword arguments to configure the synthesis.
            See :class:`SynthesisConfiguration` for a list of available options.

    Returns:
        A tuple containing the optimized circuit and the synthesis results.
    """
    config = _config_from_kwargs(kwargs)

    qc = _import_circuit(circuit)
    if initial_tableau is not None:
        if coupling_map is not None:
            synthesizer = CliffordSynthesizer(
                _import_tableau(initial_tableau, include_destabilizers), qc, set(coupling_map)
            )
        else:
            synthesizer = CliffordSynthesizer(_import_tableau(initial_tableau, include_destabilizers), qc)
    elif coupling_map is not None:
        synthesizer = CliffordSynthesizer(qc, set(coupling_map), include_destabilizers)
    else:
        synthesizer = CliffordSynthesizer(qc, include_destabilizers)

    synthesizer.synthesize(config)

    results = synthesizer.results
    circ = _circuit_from_qasm(results.circuit)

    return circ, results
