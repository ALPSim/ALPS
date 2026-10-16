# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_interop.py. Changes: pyalps names;
# spectra and reference energies come from pyalps.ed and the fock_matrix
# fixture; integral extraction is also checked on a model that uses derived
# operator names and same-site products.
"""Adapters between qdk-chemistry and ALPS representations.

The round trips must be exact: an adapter that is merely close is a silent
source of wrong physics downstream.
"""

from __future__ import annotations

import numpy as np
import pytest

from qdk_chemistry_alps import interop
from qdk_chemistry_alps.abinitio import build_abinitio_hamiltonian


def _terms(qubit_op) -> dict[str, complex]:
    d = qubit_op.to_dict()
    return {
        s: complex(re, im)
        for s, re, im in zip(
            d["pauli_strings"],
            d["coefficients"]["real"],
            d["coefficients"]["imag"],
            strict=True,
        )
    }


# --- spin models -----------------------------------------------------------

def test_qubit_operator_round_trip_is_exact():
    from qdk_chemistry.data import QubitOperator

    strings = ["IIII", "XYZI", "ZZII", "IXXI", "YIIY", "IIZZ"]
    rng = np.random.default_rng(11)
    coeffs = rng.standard_normal(len(strings))
    qo = QubitOperator(strings, coeffs)

    back = interop.hamiltonian_to_qubit_operator(
        interop.qubit_operator_to_hamiltonian(qo)
    )
    original, returned = _terms(qo), _terms(back)
    assert set(original) == set(returned)
    worst = max(abs(original[k] - returned[k]) for k in original)
    assert worst < 1e-12


def test_heisenberg_from_qdk_matches_the_dense_spectrum(fock_matrix):
    from pyalps import ed
    from qdk_chemistry.data import LatticeGraph
    from qdk_chemistry.utils.model_hamiltonians import create_heisenberg_hamiltonian

    qo = create_heisenberg_hamiltonian(
        LatticeGraph.chain(6), jx=1.0, jy=1.0, jz=1.0, hz=0.3
    )
    H = interop.qubit_operator_to_hamiltonian(qo)
    e_alps = float(ed.prepare_product(H).eigenstates(count=1).eigenvalues[0])
    e_dense = float(np.min(np.linalg.eigvalsh(fock_matrix(H))))
    assert abs(e_alps - e_dense) < 1e-10


def test_qubit_operator_matches_its_own_dense_matrix(fock_matrix):
    """Pauli characters, little-endian positions and site order, explicitly."""
    from qdk_chemistry.data import QubitOperator

    paulis = {
        "I": np.eye(2), "X": np.array([[0, 1], [1, 0]]),
        "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]),
    }
    strings, coeffs = ["XYZ", "IZX", "YII"], np.array([0.3, -0.7, 1.1])
    expected = np.zeros((8, 8), dtype=complex)
    for string, c in zip(strings, coeffs, strict=True):
        # rightmost character = qubit 0 = site 0 = most significant factor
        m = np.eye(1)
        for ch in reversed(string):
            m = np.kron(m, paulis[ch])
        expected += c * m
    H = interop.qubit_operator_to_hamiltonian(QubitOperator(strings, coeffs))
    np.testing.assert_allclose(fock_matrix(H), expected, atol=1e-14)


def test_spin_export_rejects_fermionic_sites(hubbard_alps):
    with pytest.raises(ValueError, match="dim-2 bosonic"):
        interop.hamiltonian_to_qubit_operator(hubbard_alps(2))


# --- lattices --------------------------------------------------------------

def test_lattice_graph_round_trip_preserves_adjacency():
    from qdk_chemistry.data import LatticeGraph

    n = 6
    adjacency = np.zeros((n, n))
    for i in range(n):
        j = (i + 1) % n
        w = 1.0 if i % 2 == 0 else 0.5
        adjacency[i, j] = adjacency[j, i] = w

    graph = LatticeGraph.from_dense_matrix(adjacency)
    lattice, type_weights = interop.lattice_graph_to_alps(graph)
    back = interop.alps_lattice_to_lattice_graph(lattice, type_weights)

    assert lattice.num_sites == n
    assert type_weights == [0.5, 1.0]
    assert np.allclose(np.asarray(back.adjacency_matrix()), adjacency)
    assert np.allclose(np.asarray(lattice.coordinates).ravel(), np.arange(n))


def test_graph_lattice_carries_a_model(hubbard_qdk):
    """A pyalps model on the converted lattice reproduces the QDK model."""
    from pyalps import ed, models
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import LatticeGraph

    lattice, _ = interop.lattice_graph_to_alps(LatticeGraph.chain(4))
    H = models.hubbard(lattice, t=1.0, U=4.0)
    sector = {"N": 4.0, "Sz": 0.0}
    energy = ed.prepare_product(H, sector=sector).eigenstates(count=1).eigenvalues[0]
    reference = create("multi_configuration_calculator", "macis_cas").run(
        hubbard_qdk(4, u=4.0), 2, 2
    )[0]
    assert abs(energy - reference) < 1e-10


def test_lattice_graph_rejects_asymmetric_adjacency():
    class _Asymmetric:
        def adjacency_matrix(self):
            return np.array([[0.0, 1.0], [0.0, 0.0]])

    with pytest.raises(ValueError, match="square and symmetric"):
        interop.lattice_graph_to_alps(_Asymmetric())


def test_rcm_ordering_reduces_bandwidth():
    from qdk_chemistry.data import LatticeGraph

    graph = LatticeGraph.square(3, 4)
    adjacency = np.asarray(graph.adjacency_matrix(), dtype=float)
    n = adjacency.shape[0]
    shuffled = np.random.default_rng(0).permutation(n)
    adjacency = adjacency[shuffled][:, shuffled]

    def bandwidth(a):
        i, j = np.nonzero(a)
        return int(np.max(np.abs(i - j))) if i.size else 0

    order = interop.rcm_ordering(LatticeGraph.from_dense_matrix(adjacency))
    assert sorted(order.tolist()) == list(range(n))
    assert bandwidth(adjacency[order][:, order]) < bandwidth(adjacency)


# --- fermionic integral extraction ----------------------------------------

def test_integral_extraction_inverts_the_builder(integrals):
    h1, eri = integrals(4)
    H = build_abinitio_hamiltonian(h1, eri, tol=0.0)
    h1_back, eri_back, ecore = interop.hamiltonian_to_integrals(H)

    assert np.max(np.abs(h1_back - h1)) < 1e-12
    assert np.max(np.abs(eri_back - eri)) < 1e-12
    assert abs(ecore) < 1e-12


def test_integral_extraction_reads_composite_and_derived_operators():
    """n, n_updn, Hermitian-conjugate names and a constant all decompose."""
    from pyalps import models

    site = models.spinful_fermion()

    def op(name, i):
        return models.located(site, name, i)

    H = models.Hamiltonian(models.lattice.chain(2, bc="open"), site, "strict", ["N", "Sz"])
    H.add(-0.8 * (op("c_dag_up", 0) * op("c_up", 1)), add_hc=True)
    H.add(-0.8 * (op("c_dag_dn", 0) * op("c_dn", 1)), add_hc=True)
    H.add(3.0 * op("n_updn", 0) + 2.5 * op("n_updn", 1))
    H.add(0.4 * op("n", 1) + 0.25 * op("I", 0))

    h1, eri, ecore = interop.hamiltonian_to_integrals(H)
    np.testing.assert_allclose(h1, [[0.0, -0.8], [-0.8, 0.4]], atol=1e-14)
    assert eri[0, 0, 0, 0] == pytest.approx(3.0)
    assert eri[1, 1, 1, 1] == pytest.approx(2.5)
    assert np.count_nonzero(np.abs(eri) > 1e-14) == 2
    assert ecore == pytest.approx(0.25)


@pytest.mark.parametrize(
    ("kind", "names"),
    [
        ("pairing", ("c_dag_up", "c_dag_dn")),
        ("spin-asymmetric", ("c_dag_up", "c_up")),
    ],
)
def test_integral_extraction_rejects_non_chemist_hamiltonians(kind, names):
    from pyalps import models

    site = models.spinful_fermion()
    H = models.Hamiltonian(models.lattice.chain(2, bc="open"), site, "strict")
    H.add(models.located(site, names[0], 0) * models.located(site, names[1], 1))
    with pytest.raises(ValueError):
        interop.hamiltonian_to_integrals(H)


@pytest.mark.pyscf
def test_extracted_integrals_reproduce_the_alps_energy(hubbard_alps):
    """Round trip through the integral representation preserves the spectrum."""
    from pyalps import ed
    from pyscf import fci

    H = hubbard_alps(4, u=4.0)
    h1, eri, ecore = interop.hamiltonian_to_integrals(H)

    e_alps = float(
        ed.prepare_product(H, sector={"N": 4.0, "Sz": 0.0})
        .eigenstates(count=1).eigenvalues[0]
    )
    e_pyscf = fci.direct_spin1.FCI().kernel(h1, eri, 4, (2, 2))[0] + ecore
    assert abs(e_alps - e_pyscf) < 1e-10
