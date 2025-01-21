"""Main entry point for the compilation module."""

from __future__ import annotations

from typing import TYPE_CHECKING, Union

from qiskit import QuantumCircuit, QuantumRegister, qasm3
from qiskit.transpiler import Layout, TranspileLayout

if TYPE_CHECKING:
    from os import PathLike

    from qiskit.providers import Backend
    from qiskit.providers.models import BackendProperties
    from qiskit.transpiler.target import Target

    from mqt.core.ir import QuantumComputation

    from .visualization import SearchVisualizer

    CircuitInputType = Union[QuantumComputation, str, PathLike[str], QuantumCircuit]

from mqt.core import load

from .load_architecture import load_architecture
from .load_calibration import load_calibration
from .pyqmap import (
    Arch,
    Architecture,
    CommanderGrouping,
    Configuration,
    EarlyTermination,
    Encoding,
    Heuristic,
    InitialLayout,
    Layering,
    LookaheadHeuristic,
    MappingResults,
    Method,
    SwapReduction,
    map,  # noqa: A004
)

__all__ = [
    "compile",
]


def __dir__() -> list[str]:
    return __all__


def extract_initial_layout_from_qasm(qasm: str, qregs: list[QuantumRegister]) -> Layout:
    """Extract the initial layout resulting from compiling a circuit from a QASM file.

    Args:
        qasm: The QASM file to extract the initial layout from.
        qregs: The quantum registers of the circuit.

    Returns:
        The initial layout.
    """
    for line in qasm.split("\n"):
        if line.startswith("// i "):
            # strip away initial part of line
            stripped_line = line[5:]
            # split line into tokens
            tokens = stripped_line.split(" ")
            # convert tokens to integers
            int_tokens = [int(token) for token in tokens]
            # create an empty layout
            return Layout().from_intlist(int_tokens, *qregs)
    msg = "No initial layout found in QASM file."
    raise ValueError(msg)


def compile(  # noqa: A001
    circ: CircuitInputType,
    arch: str | Arch | Architecture | Backend | None,
    calibration: str | BackendProperties | Target | None = None,
    method: str | Method = "heuristic",
    heuristic: str | Heuristic = "gate_count_max_distance",
    initial_layout: str | InitialLayout = "dynamic",
    iterative_bidirectional_routing_passes: int | None = None,
    layering: str | Layering = "individual_gates",
    automatic_layer_splits_node_limit: int | None = 5000,
    early_termination: str | EarlyTermination = "none",
    early_termination_limit: int = 0,
    lookahead_heuristic: str | LookaheadHeuristic | None = "gate_count_max_distance",
    lookaheads: int = 15,
    lookahead_factor: float = 0.5,
    use_teleportation: bool = False,
    teleportation_fake: bool = False,
    teleportation_seed: int = 0,
    encoding: str | Encoding = "commander",
    commander_grouping: str | CommanderGrouping = "fixed3",
    swap_reduction: str | SwapReduction = "coupling_limit",
    swap_limit: int = 0,
    include_WCNF: bool = False,  # noqa: N803
    use_subsets: bool = True,
    subgraph: set[int] | None = None,
    pre_mapping_optimizations: bool = True,
    post_mapping_optimizations: bool = True,
    add_measurements_to_mapped_circuit: bool = True,
    add_barriers_between_layers: bool = False,
    verbose: bool = False,
    debug: bool = False,
    visualizer: SearchVisualizer | None = None,
) -> tuple[QuantumCircuit, MappingResults]:
    """Interface to the MQT QMAP tool for mapping quantum circuits.

    Args:
        circ: The circuit to map.
        arch: The architecture to map to.
        calibration: The calibration to use.
        method: The mapping method to use. Either "heuristic" or "exact". Defaults to "heuristic".
        heuristic: The heuristic function to use for the routing search. Defaults to "gate_count_max_distance".
        initial_layout: The initial layout to use. Defaults to "dynamic".
        iterative_bidirectional_routing_passes: Number of iterative bidirectional routing passes to perform or None to disable. Defaults to None.
        layering: The layering strategy to use. Defaults to "individual_gates".
        automatic_layer_splits_node_limit: The number of expanded nodes after which to split a layer or None to disable automatic layer splitting. Defaults to 5000.
        early_termination: The early termination strategy to use, i.e. terminating the search after a goal node has been found, but before it is guarantueed to be optimal. Defaults to "none".
        early_termination_limit: The number of nodes (counted according to the early termination strategy) after which to terminate the search early. Defaults to 0.
        lookahead_heuristic: The heuristic function to use as a lookahead penalty during search or None to disable lookahead. Defaults to "gate_count_max_distance".
        lookaheads: The number of lookaheads to be used or None if no lookahead should be used. Defaults to 15.
        lookahead_factor: The rate at which the contribution of future layers to the lookahead decreases. Defaults to 0.5.
        encoding: The encoding to use for the AMO and exactly one constraints. Defaults to "naive".
        commander_grouping: The grouping strategy to use for the commander and bimander encoding. Defaults to "halves".
        swap_reduction: The swap reduction strategy to use. Defaults to "coupling_limit".
        swap_limit: Set a custom limit for max swaps per layer, for the increasing reduction strategy it sets the max swaps per layer. Defaults to 0.
        include_WCNF: Include WCNF file in the results. Defaults to False.
        use_subsets: Use qubit subsets, or consider all available physical qubits at once. Defaults to True.
        subgraph: List of qubits to consider for mapping (in exact mapper), if None all qubits are considered. Defaults to None.
        use_teleportation: Use teleportation in addition to swaps. Defaults to False.
        teleportation_fake: Assign qubits as ancillary for teleportation in the initial placement but don't actually use them (used for comparisons). Defaults to False.
        teleportation_seed: Fix a seed for the RNG in the initial ancilla placement (0 means the RNG will be seeded from /dev/urandom/ or similar). Defaults to 0.
        pre_mapping_optimizations: Run pre-mapping optimizations. Defaults to True.
        post_mapping_optimizations: Run post-mapping optimizations. Defaults to True.
        add_measurements_to_mapped_circuit: Whether to add measurements at the end of the mapped circuit. Defaults to True.
        add_barriers_between_layers: Whether to add barriers between layers to make them apparent after mapping. Defaults to False.
        verbose: Print more detailed information during the mapping process. Defaults to False.
        debug: Gather additional information during the mapping process (e.g. number of generated nodes, branching factors, ...). Defaults to False.
        visualizer: A SearchVisualizer object to log the search process to. Defaults to None.

    Returns:
        The mapped circuit and the mapping results.
    """
    if subgraph is None:
        subgraph = set()

    if arch is None and calibration is None:
        msg = "Either arch or calibration must be specified"
        raise ValueError(msg)

    architecture = load_architecture(arch)
    load_calibration(architecture, calibration)

    config = Configuration()
    config.method = Method(method)
    config.heuristic = Heuristic(heuristic)
    config.initial_layout = InitialLayout(initial_layout)
    if iterative_bidirectional_routing_passes is None:
        config.iterative_bidirectional_routing = False
    else:
        config.iterative_bidirectional_routing = True
        config.iterative_bidirectional_routing_passes = iterative_bidirectional_routing_passes
    config.layering = Layering(layering)
    if automatic_layer_splits_node_limit is None:
        config.automatic_layer_splits = False
    else:
        config.automatic_layer_splits = True
        config.automatic_layer_splits_node_limit = automatic_layer_splits_node_limit
    config.early_termination = EarlyTermination(early_termination)
    config.early_termination_limit = early_termination_limit
    config.encoding = Encoding(encoding)
    config.commander_grouping = CommanderGrouping(commander_grouping)
    config.swap_reduction = SwapReduction(swap_reduction)
    config.swap_limit = swap_limit
    config.include_WCNF = include_WCNF
    config.use_subsets = use_subsets
    config.subgraph = subgraph
    config.use_teleportation = use_teleportation
    config.teleportation_fake = teleportation_fake
    config.teleportation_seed = teleportation_seed
    config.pre_mapping_optimizations = pre_mapping_optimizations
    config.post_mapping_optimizations = post_mapping_optimizations
    config.add_measurements_to_mapped_circuit = add_measurements_to_mapped_circuit
    config.add_barriers_between_layers = add_barriers_between_layers
    config.verbose = verbose
    config.debug = debug
    if visualizer is not None and visualizer.data_logging_path is not None:
        config.data_logging_path = visualizer.data_logging_path
    if lookahead_heuristic is None:
        config.lookahead_heuristic = LookaheadHeuristic.none
        config.lookaheads = 0
    else:
        config.lookahead_heuristic = LookaheadHeuristic(lookahead_heuristic)
        config.lookaheads = lookaheads
    config.lookahead_factor = lookahead_factor

    qc = load(circ)
    results = map(qc, architecture, config)

    circ = qasm3.loads(results.mapped_circuit)
    layout = extract_initial_layout_from_qasm(results.mapped_circuit, circ.qregs)

    circ._layout = TranspileLayout(  # noqa: SLF001
        initial_layout=layout, input_qubit_mapping=layout.get_virtual_bits()
    )

    return circ, results
