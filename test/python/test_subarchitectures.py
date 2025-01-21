"""Test subarchitecture generation."""

from __future__ import annotations

import contextlib
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pytest
import rustworkx as rx
from qiskit.providers.fake_provider import GenericBackendV2

from mqt.qmap import Architecture
from mqt.qmap.subarchitectures import (
    SubarchitectureOrder,
    ibm_guadalupe_subarchitectures,
    rigetti_16_subarchitectures,
)

if TYPE_CHECKING:
    from mqt.qmap._compat.typing import TypeAlias

with contextlib.suppress(TypeError):
    Graph: TypeAlias = rx.PyGraph[int, Optional[int]]


@pytest.fixture
def ibm_guadalupe() -> SubarchitectureOrder:
    """Return the SubarchitectureOrder for the IBM Guadalupe architecture."""
    return SubarchitectureOrder.from_coupling_map([
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 5),
        (1, 4),
        (5, 8),
        (4, 7),
        (6, 7),
        (8, 9),
        (7, 10),
        (8, 11),
        (10, 12),
        (12, 15),
        (12, 13),
        (13, 14),
        (11, 14),
    ])


@pytest.fixture
def rigetti16() -> SubarchitectureOrder:
    """Return the SubarchitectureOrder for the Rigetti 16Q architecture."""
    return SubarchitectureOrder.from_coupling_map([
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 12),
        (12, 13),
        (13, 14),
        (14, 15),
        (0, 15),
        (3, 12),
        (4, 11),
    ])


@pytest.fixture
def rigetti16_opt() -> Graph:
    """Return the optimal subarchitecture candidate for the Rigetti 16Q architecture."""
    cm = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 12),
        (12, 13),
        (3, 12),
        (4, 11),
    ]

    num_nodes = max(max(int(u), int(v)) for u, v in cm)
    graph: Graph = rx.PyGraph()
    graph.add_nodes_from(list(range(num_nodes + 1)))
    graph.add_edges_from_no_data(list(cm))
    return graph


@pytest.fixture
def singleton_graph() -> Graph:
    """Return a graph with a single node."""
    g: Graph = rx.PyGraph()
    g.add_node(0)
    return g


def test_singleton_graph(singleton_graph: Graph) -> None:
    """Verify that singleton graph has trivial ordering."""
    order = SubarchitectureOrder.from_retworkx_graph(singleton_graph)

    assert len(order.sgs) == 2
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(1)[0], singleton_graph)


def test_two_node_graph(singleton_graph: Graph) -> None:
    """Verify ordering for graph with two nodes and one edge."""
    order = SubarchitectureOrder.from_coupling_map([(0, 1)])
    assert len(order.sgs) == 3
    assert len(order.sgs[0]) == 0
    assert len(order.sgs[1]) == 1
    assert len(order.sgs[2]) == 1
    assert rx.is_isomorphic(order.optimal_candidates(2)[0], order.sgs[2][0])
    assert rx.is_isomorphic(order.optimal_candidates(1)[0], singleton_graph)


def test_ibm_guadalupe_opt(ibm_guadalupe: SubarchitectureOrder) -> None:
    """Verify optimal candidates for IBM Guadalupe architecture."""
    opt_cand_9 = ibm_guadalupe.optimal_candidates(9)
    assert len(opt_cand_9) == 2
    assert opt_cand_9[0].num_nodes() == 15
    assert opt_cand_9[1].num_nodes() == 15
    assert not rx.is_isomorphic(opt_cand_9[0], opt_cand_9[1])


def test_ibm_guadalupe_cov(ibm_guadalupe: SubarchitectureOrder) -> None:
    """Verify covering for IBM Guadalupe architecture."""
    cov = ibm_guadalupe.covering(9, 2)
    assert 1 <= len(cov) <= 2

    for sg in ibm_guadalupe.sgs[9]:
        covered = False
        for co in cov:
            if rx.is_subgraph_isomorphic(co, sg):
                covered = True
                break
        assert covered


def test_rigetti16_opt(rigetti16: SubarchitectureOrder, rigetti16_opt: Graph) -> None:
    """Verify optimal candidates for Rigetti 16Q architecture."""
    opt = rigetti16.optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_rigetti16_opt_library(rigetti16_opt: Graph) -> None:
    """Verify optimal candidates for Rigetti 16Q architecture from library."""
    opt = rigetti_16_subarchitectures().optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_rigetti16_opt_library_from_str(rigetti16_opt: Graph) -> None:
    """Verify optimal candidates for Rigetti 16Q architecture from string."""
    opt = SubarchitectureOrder.from_string("rigetti_16").optimal_candidates(10)
    assert len(opt) == 1

    opt_cand = opt[0]
    assert rx.is_isomorphic(opt_cand, rigetti16_opt)


def test_ibm_guadalupe_library() -> None:
    """Verify optimal candidates for IBM Guadalupe architecture from library."""
    opt_cand_9 = ibm_guadalupe_subarchitectures().optimal_candidates(9)
    assert len(opt_cand_9) == 2
    assert opt_cand_9[0].num_nodes() == 15
    assert opt_cand_9[1].num_nodes() == 15
    assert not rx.is_isomorphic(opt_cand_9[0], opt_cand_9[1])


def test_store_subarch(ibm_guadalupe: SubarchitectureOrder) -> None:
    """Verify that subarchitecture order can be stored and loaded."""
    ibm_guadalupe.store_library("tmp")

    p = Path("tmp.pickle")

    loaded_tmp = SubarchitectureOrder.from_library(p)

    if p.exists():
        p.unlink()

    opt_origin = ibm_guadalupe.optimal_candidates(8)
    opt_loaded = loaded_tmp.optimal_candidates(8)

    assert len(opt_origin) == len(opt_loaded)
    for opt_cand_orig, opt_cand_load in zip(opt_origin, opt_loaded):
        assert rx.is_isomorphic(opt_cand_load, opt_cand_orig)


def test_subarchitecture_from_qmap_arch() -> None:
    """Verify that subarchitecture order can be created from QMAP architectures."""
    cm = {(0, 1), (1, 0), (1, 2), (2, 1)}
    arch = Architecture(3, cm)
    so_arch = SubarchitectureOrder.from_qmap_architecture(arch)
    so_cm = SubarchitectureOrder.from_coupling_map(cm)

    assert so_arch.subarch_order == so_cm.subarch_order


def test_subarchitecture_from_qiskit_backend() -> None:
    """Verify that subarchitecture order can be created from Qiskit backends."""
    arch = GenericBackendV2(num_qubits=5, coupling_map=[[0, 1], [1, 0], [1, 2], [2, 1], [1, 3], [3, 1], [3, 4], [4, 3]])
    so_arch = SubarchitectureOrder.from_backend_v2(arch)
    cm = [(c[0], c[1]) for c in arch.coupling_map]
    so_cm = SubarchitectureOrder.from_coupling_map(cm)

    assert so_arch.subarch_order == so_cm.subarch_order


def test_invalid_opt_cand_arg(ibm_guadalupe: SubarchitectureOrder) -> None:
    """Verify that invalid arguments for optimal candidates raise an error."""
    with pytest.raises(
        ValueError,
        match=r"Number of qubits must not be smaller or equal 0 or larger then number of physical qubits of architecture.",
    ):
        ibm_guadalupe.optimal_candidates(100)
