# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# WIP/classic-dmrg-chemistry at d4384a7f8,
# python/qdk-chemistry-alps/tests/test_classic_calculator.py. Changes: alpsx
# compared the classic RDMs with its MPS calculator, which is not ported;
# here they are compared element by element with PySCF's FCI RDMs, and the
# energies with alps_ed and MACIS. Written here: the CI-vector, container
# selection and measurement tests.
"""The classic DMRG calculator: energies, RDM conventions, the CI vector,
the containers and the measurements.

Tests that need ``RdmContainer`` are marked ``local_qdk``; the rest run on
released qdk-chemistry, where the calculator returns a state vector.
"""

from __future__ import annotations

from math import comb

import numpy as np
import pytest
from qdk_chemistry.algorithms import create

import qdk_chemistry_alps  # noqa: F401

pytestmark = [pytest.mark.classic_dmrg]

CONVERGED = {"max_bond_dim": 128, "warmup_states": 128, "sweeps": 5,
             "convergence_tol": 1e-12}


def _classic(**settings):
    return create("multi_configuration_calculator", "alps_classic_dmrg", **settings)


def _dense_hamiltonian(h1, eri, ecore):
    from qdk_chemistry.data import (
        CanonicalFourCenterHamiltonianContainer,
        Hamiltonian,
        ModelOrbitals,
    )

    n = h1.shape[0]
    return Hamiltonian(CanonicalFourCenterHamiltonianContainer(
        h1, eri.reshape(-1), ModelOrbitals(n), ecore, h1))


def _assert_rdms_match_fci(calculator, wavefunction, h1, eri, n_alpha, n_beta):
    from pyscf import fci

    n = h1.shape[0]
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-14
    _, ci = solver.kernel(h1, eri, n, (n_alpha, n_beta))
    dm1, dm2 = solver.make_rdm12(ci, n, (n_alpha, n_beta))
    (dm1a, dm1b), (dm2aa, dm2ab, dm2bb) = solver.make_rdm12s(ci, n, (n_alpha, n_beta))
    got1, got2 = calculator.make_rdms()
    np.testing.assert_allclose(got1, dm1, atol=1e-7)
    np.testing.assert_allclose(got2, dm2, atol=1e-7)
    np.testing.assert_allclose(calculator.last_result.one_rdm_alpha, dm1a, atol=1e-7)
    np.testing.assert_allclose(calculator.last_result.one_rdm_beta, dm1b, atol=1e-7)
    # The container's spin blocks, in QDK's aaaa/aabb/bbbb layout.
    blocks = [np.asarray(b).reshape((n,) * 4) for b in wavefunction.get_active_two_rdm_spin_dependent()]
    for got, expected in zip(blocks, (dm2aa, dm2ab, dm2bb), strict=True):
        np.testing.assert_allclose(got, expected, atol=1e-7)
    np.testing.assert_allclose(np.asarray(wavefunction.get_active_two_rdm_spin_traced()).reshape((n,) * 4),
                               dm2, atol=1e-7)


def _coefficients(wavefunction) -> dict[str, float]:
    return {det.to_string(): float(c) for det, c in zip(
        wavefunction.get_active_determinants(), wavefunction.get_coefficients(), strict=True)}


def _assert_same_state(ours, reference, atol=1e-8):
    """Equal coefficients determinant by determinant, up to one global sign."""
    significant = {det for det, c in reference.items() if abs(c) > 1e-10}
    assert significant <= ours.keys()
    sign = np.sign(sum(reference[det] * ours[det] for det in significant))
    for det in ours.keys() | reference.keys():
        assert sign * ours.get(det, 0.0) == pytest.approx(reference.get(det, 0.0), abs=atol), det


def _sector_problem(h1, eri, n_alpha, n_beta):
    """pyalps.ed's ground state of the integrals and its sector problem."""
    from qdk_chemistry_alps.abinitio import build_abinitio_hamiltonian, solve_ground_state

    H = build_abinitio_hamiltonian(h1, eri)
    return H, solve_ground_state(H, n_alpha, n_beta)


def _alps_vector(wavefunction, states):
    """A QDK state vector on the rows of a pyalps.ed sector problem."""
    from qdk_chemistry_alps.abinitio import interleaved_to_blocked_parity

    local = {"0": 0, "u": 1, "d": 2, "2": 3}
    states = np.asarray(states)
    row = {tuple(int(x) for x in state): i for i, state in enumerate(states)}
    psi = np.zeros(len(row))
    for det, c in _coefficients(wavefunction).items():
        occupations = [local[x] for x in det[:states.shape[1]]]  # to_string pads to the capacity
        alpha = np.array([x & 1 for x in occupations])
        beta = np.array([x >> 1 for x in occupations])
        psi[row[tuple(occupations)]] = c * interleaved_to_blocked_parity(alpha, beta)
    return psi


@pytest.mark.local_qdk
@pytest.mark.parametrize(("n", "na", "nb"), [(2, 1, 1), (3, 2, 1), (4, 2, 2)])
def test_classic_matches_ed_and_macis(hubbard_qdk, n, na, nb):
    h = hubbard_qdk(n, u=2.0)
    classic = _classic(**CONVERGED, truncation_tol=1e-14, wavefunction_container="rdm")
    e_classic, wavefunction = classic.run(h, na, nb)
    e_ed, _ = create("multi_configuration_calculator", "alps_ed").run(h, na, nb)
    e_cas, _ = create("multi_configuration_calculator", "macis_cas").run(h, na, nb)
    assert e_classic == pytest.approx(e_ed, abs=1e-8)
    assert e_classic == pytest.approx(e_cas, abs=1e-8)
    assert wavefunction.get_container_type() == "rdm"
    stored = wavefunction.get_container().get_active_one_rdm_spin_traced()
    np.testing.assert_allclose(stored, classic.make_rdms()[0], atol=1e-12)
    np.testing.assert_allclose(classic.last_energies, [e_classic])
    one, two = classic.make_rdms()
    assert np.trace(one) == pytest.approx(na + nb, abs=1e-8)
    assert np.einsum("ppqq->", two) == pytest.approx((na + nb) * (na + nb - 1), abs=1e-7)


@pytest.mark.pyscf
@pytest.mark.parametrize(("n", "na", "nb"), [(3, 2, 1), (4, 2, 2)])
def test_classic_rdms_match_pyscf_fci_element_by_element(integrals, n, na, nb):
    """Index conventions of dm1 and dm2, which an energy contraction cannot see."""
    h1, eri = integrals(n)
    calculator = _classic(**CONVERGED)
    energy, wavefunction = calculator.run(_dense_hamiltonian(h1, eri, 0.713), na, nb)
    assert wavefunction.get_container_type() == "state_vector"
    from pyscf import fci

    reference = fci.direct_spin1.FCI().kernel(h1, eri, n, (na, nb))[0] + 0.713
    assert energy == pytest.approx(reference, abs=1e-9)
    _assert_rdms_match_fci(calculator, wavefunction, h1, eri, na, nb)


@pytest.mark.local_qdk
@pytest.mark.pyscf
def test_classic_rdm_container_carries_the_spin_blocks(integrals):
    h1, eri = integrals(4)
    calculator = _classic(**CONVERGED, wavefunction_container="rdm")
    _, wavefunction = calculator.run(_dense_hamiltonian(h1, eri, 0.0), 2, 1)
    assert wavefunction.get_container_type() == "rdm"
    assert calculator.last_result.ci_coefficients.size == 0
    _assert_rdms_match_fci(calculator, wavefunction, h1, eri, 2, 1)


def test_classic_requires_successful_run_for_rdms():
    for read in (lambda c: c.make_rdms(), lambda c: c.measure_local("n"),
                 lambda c: c.measure_correlations("Sz"), lambda c: c.spin_spin_correlations()):
        with pytest.raises(RuntimeError, match="has not produced"):
            read(_classic())


@pytest.mark.parametrize(("n", "na", "nb"), [(2, 1, 1), (3, 2, 1), (4, 2, 2), (5, 3, 2)])
def test_classic_ci_vector_matches_ed_and_macis(hubbard_qdk, n, na, nb):
    """Determinant labels and signs, against both solvers, up to one global sign."""
    h = hubbard_qdk(n, u=2.0)
    classic = _classic(**CONVERGED)
    energy, wavefunction = classic.run(h, na, nb)
    assert wavefunction.get_container_type() == "state_vector"
    ours = _coefficients(wavefunction)
    assert len(ours) == comb(n, na) * comb(n, nb)
    assert sum(c * c for c in ours.values()) == pytest.approx(1.0, abs=1e-12)
    for name in ("alps_ed", "macis_cas"):
        reference_energy, reference = create("multi_configuration_calculator", name).run(h, na, nb)
        assert energy == pytest.approx(reference_energy, abs=1e-8)
        _assert_same_state(ours, _coefficients(reference))
    # The run's RDMs travel with the vector.
    one, two = classic.make_rdms()
    np.testing.assert_allclose(wavefunction.get_active_one_rdm_spin_traced(), one, atol=1e-12)
    alpha, beta = wavefunction.get_active_one_rdm_spin_dependent()
    np.testing.assert_allclose(alpha, classic.last_result.one_rdm_alpha, atol=1e-12)
    np.testing.assert_allclose(beta, classic.last_result.one_rdm_beta, atol=1e-12)
    np.testing.assert_allclose(np.asarray(wavefunction.get_active_two_rdm_spin_traced()).reshape((n,) * 4),
                               two, atol=1e-12)


def test_truncated_ci_vector_is_the_dmrg_state(integrals):
    """Far from FCI, <psi|H|psi> of the returned vector is the DMRG energy."""
    h1, eri = integrals(6, seed=8)
    calculator = _classic(max_bond_dim=6, warmup_states=6, sweeps=2, convergence_tol=1e-12)
    energy, wavefunction = calculator.run(_dense_hamiltonian(h1, eri, 0.25), 3, 2)
    H, exact = _sector_problem(h1, eri, 3, 2)
    assert energy - 0.25 > exact.eigenvalues[0] + 1e-4  # truncation matters here
    psi = _alps_vector(wavefunction, exact.problem.states)
    matrix = exact.problem.operator_matrix(H.operator_decomposition().op).to_scipy()
    assert psi @ psi == pytest.approx(1.0, abs=1e-12)
    assert psi @ (matrix @ psi) + 0.25 == pytest.approx(energy, abs=1e-9)
    # The run's RDMs are those of the same state.
    dm1, dm2 = calculator.make_rdms()
    assert np.einsum("pq,pq->", h1, dm1) + 0.5 * np.einsum("pqrs,pqrs->", eri, dm2) + 0.25 == pytest.approx(
        energy, abs=1e-9)


def test_spin_blocks_follow_qdks_layout(integrals):
    """QDK's aaaa/aabb/bbbb blocks, as MACIS fills them, are two_rdm_aa/ab/bb."""
    h1, eri = integrals(3)
    h = _dense_hamiltonian(h1, eri, 0.0)
    calculator = _classic(**CONVERGED)
    calculator.run(h, 2, 1)
    macis = create("multi_configuration_calculator", "macis_cas", calculate_one_rdm=True, calculate_two_rdm=True)
    _, reference = macis.run(h, 2, 1)
    aaaa, aabb, bbbb = (np.asarray(b).reshape((3,) * 4) for b in reference.get_active_two_rdm_spin_dependent())
    result = calculator.last_result
    np.testing.assert_allclose(result.two_rdm_aa, aaaa, atol=1e-8)
    np.testing.assert_allclose(result.two_rdm_ab, aabb, atol=1e-8)
    np.testing.assert_allclose(result.two_rdm_bb, bbbb, atol=1e-8)
    # The mixed block is not symmetric under the swap of its electron pairs.
    assert np.abs(result.two_rdm_ab.transpose(2, 3, 0, 1) - aabb).max() > 1e-2


def test_measurements_match_exact_diagonalization(integrals):
    from pyalps import models

    h1, eri = integrals(4, seed=6)
    calculator = _classic(**CONVERGED)
    fs = models.spinful_fermion()
    at = lambda op, i: models.located(fs, op, i)  # noqa: E731
    calculator.observables = {"pair hop": at("c_dag_up", 0) * at("c_dag_dn", 0) * at("c_dn", 3) * at("c_up", 3)
                              + at("c_dag_up", 3) * at("c_dag_dn", 3) * at("c_dn", 0) * at("c_up", 0)}
    _, wavefunction = calculator.run(_dense_hamiltonian(h1, eri, 0.0), 2, 2)
    H, exact = _sector_problem(h1, eri, 2, 2)
    psi = np.asarray(exact.eigenvectors)[:, 0]

    def expect(term):
        return psi @ (exact.problem.operator_matrix(H.operator_decomposition(term).op).to_scipy() @ psi)

    n = 4
    for op in ("n", "Sz", "n_updn"):
        np.testing.assert_allclose(calculator.measure_local(op), [expect(at(op, i)) for i in range(n)], atol=1e-8)
    for op_i, op_j, left, right in (("Sz", None, "Sz", "Sz"), ("n", "Sz", "n", "Sz"),
                                    ("Sp", None, "Sm", "Sp"), ("c_dn", None, "c_dag_dn", "c_dn")):
        expected = [[expect(at(left, i) * at(right, j)) for j in range(n)] for i in range(n)]
        np.testing.assert_allclose(calculator.measure_correlations(op_i, op_j), expected, atol=1e-8)
    spin = calculator.spin_spin_correlations()
    assert spin.sum() == pytest.approx(wavefunction.get_container().compute_s_squared(), abs=1e-8)
    assert calculator.last_result.observables["pair hop"] == pytest.approx(
        expect(calculator.observables["pair hop"]), abs=1e-8)
    with pytest.raises(ValueError, match="calculator.observables"):
        calculator.measure_correlations("c_up", "n")


def test_ci_cutoff_trims_the_vector_and_drops_the_rdms(hubbard_qdk):
    h = hubbard_qdk(4, u=4.0)
    full = _classic(**CONVERGED).run(h, 2, 2)[1]
    trimmed = _classic(**CONVERGED, ci_coefficient_cutoff=0.05).run(h, 2, 2)[1]
    assert len(trimmed.get_coefficients()) < len(full.get_coefficients())
    assert np.linalg.norm(trimmed.get_coefficients()) == pytest.approx(1.0, abs=1e-12)
    assert full.get_container().has_one_rdm_spin_traced()
    assert not trimmed.get_container().has_one_rdm_spin_traced()


def test_container_choice_is_made_before_the_solve(hubbard_qdk):
    h = hubbard_qdk(4)
    calculator = _classic(wavefunction_container="statevector", ci_memory_limit_gb=1e-9)
    with pytest.raises(MemoryError, match="ci_memory_limit_gb"):
        calculator.run(h, 2, 2)
    assert calculator.last_result is None
    with pytest.raises(ValueError, match="wavefunction_container"):
        _classic(wavefunction_container="invalid").run(h, 2, 2)


@pytest.mark.local_qdk
def test_auto_falls_back_to_the_rdm_container(hubbard_qdk):
    calculator = _classic(sweeps=4, ci_memory_limit_gb=1e-9)
    _, wavefunction = calculator.run(hubbard_qdk(4), 2, 2)
    assert wavefunction.get_container_type() == "rdm"
    assert calculator.last_result.ci_coefficients.size == 0


@pytest.mark.local_qdk
def test_classic_rdm_container_roundtrip(hubbard_qdk, tmp_path):
    from qdk_chemistry.data import Wavefunction

    calculator = _classic(sweeps=4, wavefunction_container="rdm")
    _, wavefunction = calculator.run(hubbard_qdk(2), 1, 1)
    path = str(tmp_path / "classic.wavefunction.json")
    wavefunction.to_json_file(path)
    restored = Wavefunction.from_json_file(path)
    assert restored.get_container_type() == "rdm"
    assert restored.get_container().get_metadata() == wavefunction.get_container().get_metadata()
    np.testing.assert_allclose(restored.get_container().get_active_one_rdm_spin_traced(),
                               calculator.make_rdms()[0], atol=1e-12)


def test_dense_molecular_integrals_match_ed(integrals):
    h1, eri = integrals(4)
    h = _dense_hamiltonian(h1, eri, .713)
    classic = _classic(max_bond_dim=128, sweeps=8, convergence_tol=1e-12,
                       truncation_tol=1e-14, warmup_states=128)
    e_classic, _ = classic.run(h, 2, 1)
    e_ed, _ = create("multi_configuration_calculator", "alps_ed").run(h, 2, 1)
    assert e_classic == pytest.approx(e_ed, abs=1e-9)


def test_unrestricted_hamiltonian_is_rejected():
    class Unrestricted:
        def is_restricted(self):
            return False

    with pytest.raises(NotImplementedError, match="restricted"):
        _classic()._run_impl(Unrestricted(), 1, 1)


@pytest.mark.pyscf
def test_lih_default_fixed_state_schedule_matches_fci():
    from pyscf import ao2mo, fci, gto, scf

    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g", verbose=0)
    mean_field = scf.RHF(mol)
    mean_field.conv_tol = 1e-12
    mean_field.kernel()
    coefficients = mean_field.mo_coeff
    n = coefficients.shape[1]
    h1 = coefficients.T @ mean_field.get_hcore() @ coefficients
    eri = ao2mo.restore(1, ao2mo.kernel(mol, coefficients), n)
    reference = fci.direct_spin1.FCI()
    e_fci, ci = reference.kernel(h1, eri, n, (2, 2), ecore=mol.energy_nuc())
    rdms = reference.make_rdm12(ci, n, (2, 2))
    h = _dense_hamiltonian(h1, eri, mol.energy_nuc())
    calculator = _classic(max_bond_dim=128, sweeps=6, convergence_tol=1e-12)
    assert calculator.settings().get("truncation_tol") == 0.0
    energy, _ = calculator.run(h, 2, 2)
    assert energy == pytest.approx(e_fci, abs=1e-9)
    for got, expected in zip(calculator.make_rdms(), rdms, strict=True):
        np.testing.assert_allclose(got, expected, atol=1e-7)
