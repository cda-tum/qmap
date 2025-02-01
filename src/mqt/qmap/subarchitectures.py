"""Functionality for computing good subarchitectures for quantum circuit mapping.

This file implements the methods presented in https://arxiv.org/abs/2210.09321.
"""

from __future__ import annotations

import pickle
from itertools import combinations
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from ._compat.importlib import resources

if TYPE_CHECKING:
    from collections.abc import Iterable

    from matplotlib import figure
    from qiskit.providers import BackendV1, BackendV2

    from . import Architecture
    from ._compat.typing import TypeAlias

import contextlib

import rustworkx as rx
import rustworkx.visualization as rxviz

__all__ = [
    "SubarchitectureOrder",
    "ibm_guadalupe_subarchitectures",
    "rigetti_16_subarchitectures",
]


def __dir__() -> list[str]:
    return __all__


with contextlib.suppress(TypeError):
    Graph: TypeAlias = rx.PyGraph[int, Optional[int]]

PartialOrder: TypeAlias = dict[tuple[int, int], set[tuple[int, int]]]

#: Architectures for which precomputed orderings are available
precomputed_backends = ["rigetti_16", "ibm_guadalupe_16"]


class SubarchitectureOrder:
    """Class representing the partial order of (sub)architectures."""

    inactive_color: str = "#1f78b4"
    active_color: str = "#faf18e"

    def __init__(self) -> None:
        """Initialize a partial order."""
        self.arch: Graph = rx.PyGraph()
        self.subarch_order: PartialOrder = {}
        self.desirable_subarchitectures: PartialOrder = {}
        self.isomorphisms: dict[tuple[int, int], dict[tuple[int, int], dict[int, int]]] = {}

        self.__compute_subarchs()
        self.__compute_subarch_order()
        self.__compute_desirable_subarchitectures()

    @classmethod
    def from_retworkx_graph(cls, graph: Graph) -> SubarchitectureOrder:
        """Construct the partial order from retworkx graph.

        Args:
            graph: retworkx graph representing the architecture.

        Returns:
            The resulting partial order.
        """
        so = SubarchitectureOrder()
        so.arch = graph

        so.__compute_subarchs()  # noqa: SLF001
        so.__compute_subarch_order()  # noqa: SLF001
        so.__compute_desirable_subarchitectures()  # noqa: SLF001
        return so

    @classmethod
    def from_coupling_map(cls, coupling_map: Iterable[tuple[int, int]]) -> SubarchitectureOrder:
        """Construct partial order from coupling map defined as set of tuples of connected qubits.

        Args:
            coupling_map: Iterable of tuples of connected qubits.

        Returns:
            The resulting partial order.
        """
        num_nodes = max(max(int(u), int(v)) for u, v in coupling_map)
        graph: Graph = rx.PyGraph()
        graph.add_nodes_from(list(range(num_nodes + 1)))
        graph.add_edges_from_no_data(list(coupling_map))

        return cls.from_retworkx_graph(graph)

    @classmethod
    def from_backend(cls, backend: BackendV1) -> SubarchitectureOrder:
        """Construct the partial order from a coupling map defined as a Qiskit backend.

        Args:
            backend: Qiskit backend.

        Returns:
            The resulting partial order.
        """
        coupling_map = [(c[0], c[1]) for c in backend.configuration().coupling_map]
        return cls.from_coupling_map(coupling_map)

    @classmethod
    def from_backend_v2(cls, backend: BackendV2) -> SubarchitectureOrder:
        """Construct the partial order from a coupling map defined as a Qiskit backend.

        Args:
            backend: Qiskit backend.

        Returns:
            The resulting partial order.
        """
        coupling_map = [(c[0], c[1]) for c in backend.coupling_map]
        return cls.from_coupling_map(coupling_map)

    @classmethod
    def from_qmap_architecture(cls, arch: Architecture) -> SubarchitectureOrder:
        """Construct the partial order from a QMAP :class:`~mqt.qmap.Architecture` object.

        Args:
            arch: QMAP architecture.

        Returns:
            The resulting partial order.
        """
        return cls.from_coupling_map(arch.coupling_map)

    @classmethod
    def from_library(cls, lib_name: str | Path) -> SubarchitectureOrder:
        """Construct the partial order from a stored library.

        Args:
            lib_name: Path to the library.

        Returns:
            The resulting partial order.
        """
        path = Path(lib_name).with_suffix(".pickle")
        with path.open("rb") as f:
            temp = pickle.load(f)  # noqa: S301

        so = SubarchitectureOrder()
        so.__dict__.update(temp.__dict__)

        return so

    @classmethod
    def from_string(cls, name: str) -> SubarchitectureOrder:
        """Construct the partial order from a library name.

        Args:
            name: Name of the library.

        Returns:
            The resulting partial order.
        """
        if name in precomputed_backends:
            ref = resources.files("mqt.qmap") / "libs" / (name + ".pickle")
            with resources.as_file(ref) as lib_path:
                return cls.from_library(lib_path)
        return SubarchitectureOrder()

    def optimal_candidates(self, nqubits: int) -> list[Graph]:
        """Return optimal subarchitecture candidate.

        Args:
            nqubits:
                size of circuit for which the optimal candidate should be given.

        Returns:
            List of optimal subarchitecture candidates for circuits of the given size.
        """
        if nqubits <= 0 or nqubits > self.arch.num_nodes():
            msg = "Number of qubits must not be smaller or equal 0 or larger then number of physical qubits of architecture."
            raise ValueError(msg)

        if nqubits == self.arch.num_nodes():
            return [self.arch]

        cands = self.__cand(nqubits)
        trans_ord = self.__transitive_closure(self.subarch_order)
        ref_ord = self.__reflexive_closure(trans_ord)

        opt_cands = set(ref_ord[next(iter(cands))])
        for cand in cands:
            opt_cands = opt_cands.intersection(set(ref_ord[cand]))

        ordered_cands = list(opt_cands)
        ordered_cands.sort()
        for cand in ordered_cands:
            opt_cands = opt_cands.difference(trans_ord[cand])

        return [self.sgs[n][i] for (n, i) in opt_cands]

    def covering(self, nqubits: int, size: int) -> list[Graph]:
        """Return covering for nqubit circuits.

        Args:
            nqubits:
                size of circuit for which the covering should be given.
            size:
                limit for the size of the covering.

        Returns:
            Subarchitecture covering for circuits of the given size.

        Note:
            A smaller covering might be found.
        """
        cov = self.__cand(nqubits)
        po_trans = self.__transitive_closure(self.subarch_order)
        ref_trans_po = self.__reflexive_closure(po_trans)
        queue = list({el for cand in cov for el in ref_trans_po[cand]})
        queue.sort(reverse=True)

        po_inv = self.__inverse_relation(po_trans)

        while len(cov) > size:
            d = queue.pop()
            cov_d = cov.intersection(po_inv[d])
            if len(cov_d) > 1:
                cov = cov.difference(cov_d)
                cov.add(d)

        return [self.sgs[n][i] for n, i in cov]

    def store_library(self, lib_name: str | Path) -> None:
        """Store ordering.

        Args:
            lib_name: Path to the library.
        """
        path = Path(lib_name).with_suffix(".pickle")
        with path.open("wb") as f:
            pickle.dump(self, f)

    def draw_subarchitecture(self, subarchitecture: Graph | tuple[int, int]) -> figure.Figure:
        """Create a matplotlib figure showing subarchitecture within the entire architecture.

        Nodes that are part of the subarchitecture are drawn yellow.
        Nodes that are not part of the subarchitecture are drawn blue.

        Args:
            subarchitecture: Subarchitecture to be drawn.

        Returns:
            Matplotlib figure.
        """
        if isinstance(subarchitecture, tuple):
            subarchitecture = self.sgs[subarchitecture[0]][subarchitecture[1]]
        colors = [SubarchitectureOrder.inactive_color for _ in range(self.arch.num_nodes())]
        for node in subarchitecture.nodes():
            colors[node] = SubarchitectureOrder.active_color
        return rxviz.mpl_draw(self.arch, node_color=colors)

    def draw_subarchitectures(self, subarchitectures: list[Graph] | list[tuple[int, int]]) -> list[figure.Figure]:
        """Create matplotlib figures showing subarchitectures within the entire architecture.

        For each subarchitecture one figure is drawn.
        Nodes that are part of the subarchitecture are drawn yellow.
        Nodes that are not part of the subarchitecture are drawn blue.

        Args:
            subarchitectures: Subarchitectures to be drawn.

        Returns:
            List of matplotlib figures.
        """
        return [self.draw_subarchitecture(subarchitecture) for subarchitecture in subarchitectures]

    def __compute_subarchs(self) -> None:
        """Compute all subarchitectures of the architecture."""
        self.sgs: list[list[Graph]] = [[] for i in range(self.arch.num_nodes() + 1)]

        for i in range(1, self.arch.num_nodes() + 1):
            node_combinations = combinations(range(self.arch.num_nodes()), i)
            for sg in (self.arch.subgraph(selected_nodes) for selected_nodes in node_combinations):
                if rx.is_connected(sg):
                    new_class = True
                    for g in self.sgs[i]:
                        if rx.is_isomorphic(g, sg):
                            new_class = False
                            break
                    if new_class:
                        self.sgs[i].append(sg)
        # init orders
        for n in range(self.arch.num_nodes() + 1):
            for i in range(len(self.sgs[n])):
                self.subarch_order[n, i] = set()
                self.desirable_subarchitectures[n, i] = set()
                self.isomorphisms[n, i] = {}

    def __compute_subarch_order(self) -> None:
        """Compute subarchitecture order."""
        for n, sgs_n in enumerate(self.sgs[:-1]):
            for i, sg in enumerate(sgs_n):
                for j, parent_sg in enumerate(self.sgs[n + 1]):
                    matcher = rx.graph_vf2_mapping(parent_sg, sg, subgraph=True)
                    for iso in matcher:
                        self.subarch_order[n, i].add((n + 1, j))
                        iso_rev = {val: key for key, val in iso.items()}
                        self.isomorphisms[n, i][n + 1, j] = iso_rev
                        break  # One isomorphism suffices

    def __complete_isos(self) -> None:
        """Complete isomorphisms."""
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                for _, i_prime in self.subarch_order[n, i]:
                    self.__combine_iso_with_parent(n, i, i_prime)

    def __combine_iso_with_parent(self, n: int, i: int, j: int) -> None:
        """Combine all isomorphisms from sgs[n][i] with those from sgs[n+1][j]."""
        first = self.isomorphisms[n, i][n + 1, j]
        for (row, k), second in self.isomorphisms[n + 1, j].items():
            self.isomorphisms[n, i][row, k] = SubarchitectureOrder.__combine_isos(first, second)

    @staticmethod
    def __combine_isos(first: dict[int, int], second: dict[int, int]) -> dict[int, int]:
        """Combine two isomorphisms."""
        combined = {}
        for src, img in first.items():
            combined[src] = second[img]
        return combined

    def __transitive_closure(self, po: PartialOrder) -> PartialOrder:
        """Compute transitive closure of partial order."""
        po_trans: PartialOrder = {(self.arch.num_nodes(), 0): set()}

        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                new_rel = set(po[n, i])
                po_trans[n, i] = new_rel.copy()
                for n_prime, i_prime in po_trans[n, i]:
                    new_rel = new_rel.union(po_trans[n_prime, i_prime])
                po_trans[n, i] = new_rel

        return po_trans

    @classmethod
    def __reflexive_closure(cls, po: PartialOrder) -> PartialOrder:
        """Compute reflexive closure of partial order."""
        po_ref = {}
        for k, v in po.items():
            v_copy = v.copy()
            v_copy.add(k)
            po_ref[k] = v_copy
        return po_ref

    def __inverse_relation(self, po: PartialOrder) -> PartialOrder:
        """Compute inverse relation of partial order."""
        po_inv: PartialOrder = {}
        for n in range(self.arch.num_nodes() + 1):
            for i in range(len(self.sgs[n])):
                po_inv[n, i] = set()
        for k, v in po.items():
            for e in v:
                po_inv[e].add(k)
        return po_inv

    def __path_order_less(self, n: int, i: int, n_prime: int, i_prime: int) -> bool:
        """Check if sgs[n][i] is less than sgs[n_prime][i_prime] in the path order."""
        lhs = self.sgs[n][i]
        rhs = self.sgs[n_prime][i_prime]
        iso = self.isomorphisms[n, i][n_prime, i_prime]
        for v in range(lhs.num_nodes()):
            for w in range(lhs.num_nodes()):
                if v is w:
                    continue
                if (
                    rx.dijkstra_shortest_path_lengths(lhs, v, lambda _x: 1, goal=w)[w]
                    > rx.dijkstra_shortest_path_lengths(rhs, iso[v], lambda _x: 1, goal=iso[w])[iso[w]]
                ):
                    return True
        return False

    def __compute_desirable_subarchitectures(self) -> None:
        """Compute desirable subarchitectures."""
        self.__complete_isos()
        for n in reversed(range(1, len(self.sgs[:-1]))):
            for i in range(len(self.sgs[n])):
                val = self.isomorphisms[n, i]
                for n_prime, i_prime in val:
                    if self.__path_order_less(n, i, n_prime, i_prime):
                        self.desirable_subarchitectures[n, i].add((n_prime, i_prime))
                des = list(self.desirable_subarchitectures[n, i])
                des.sort()
                new_des: set[tuple[int, int]] = set()
                for j, (n_prime, i_prime) in enumerate(reversed(des)):
                    idx = len(des) - j - 1
                    if not any((n_prime, i_prime) in self.subarch_order[k] for k in des[:idx]):
                        new_des.add((n_prime, i_prime))

                self.desirable_subarchitectures[n, i] = new_des
                if len(self.desirable_subarchitectures[n, i]) == 0:
                    self.desirable_subarchitectures[n, i].add((n, i))
        self.desirable_subarchitectures[self.arch.num_nodes(), 0] = {(self.arch.num_nodes(), 0)}

    def __cand(self, nqubits: int) -> set[tuple[int, int]]:
        return {
            des for (n, i), desirables in self.desirable_subarchitectures.items() if n == nqubits for des in desirables
        }


def ibm_guadalupe_subarchitectures() -> SubarchitectureOrder:
    """Load the precomputed ibm guadalupe subarchitectures.

    Returns:
        The subarchitecture order for the ibm_guadalupe architecture.
    """
    ref = resources.files("mqt.qmap") / "libs" / "ibm_guadalupe_16.pickle"
    with resources.as_file(ref) as path:
        return SubarchitectureOrder.from_library(path)


def rigetti_16_subarchitectures() -> SubarchitectureOrder:
    """Load the precomputed rigetti subarchitectures.

    Returns:
        The subarchitecture order for the 16-qubit Rigetti architecture.
    """
    ref = resources.files("mqt.qmap") / "libs" / "rigetti_16.pickle"
    with resources.as_file(ref) as path:
        return SubarchitectureOrder.from_library(path)
