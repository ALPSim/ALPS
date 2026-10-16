# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# WIP/classic-dmrg-chemistry at d4384a7f8,
# python/qdk-chemistry-alps/tests/test_classic_calculator.py. Changes: alpsx
# compared the classic RDMs with its MPS calculator, which is not ported;
# here they are compared element by element with PySCF's FCI RDMs, and the
# energies with alps_ed and MACIS.
"""The classic DMRG calculator: energies, RDM conventions and the container."""

from __future__ import annotations

import numpy as np
import pytest
from qdk_chemistry.algorithms import create

import qdk_chemistry_alps  # noqa: F401

pytestmark = [pytest.mark.classic_dmrg, pytest.mark.local_qdk]

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


def _assert_rdms_match_fci(calculator, h1, eri, n_alpha, n_beta):
    from pyscf import fci

    n = h1.shape[0]
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-14
    _, ci = solver.kernel(h1, eri, n, (n_alpha, n_beta))
    dm1, dm2 = solver.make_rdm12(ci, n, (n_alpha, n_beta))
    (dm1a, dm1b) = solver.make_rdm1s(ci, n, (n_alpha, n_beta))
    got1, got2 = calculator.make_rdms()
    np.testing.assert_allclose(got1, dm1, atol=1e-7)
    np.testing.assert_allclose(got2, dm2, atol=1e-7)
    np.testing.assert_allclose(calculator.last_result.one_rdm_alpha, dm1a, atol=1e-7)
    np.testing.assert_allclose(calculator.last_result.one_rdm_beta, dm1b, atol=1e-7)


@pytest.mark.parametrize(("n", "na", "nb"), [(2, 1, 1), (3, 2, 1), (4, 2, 2)])
def test_classic_matches_ed_and_macis(hubbard_qdk, n, na, nb):
    h = hubbard_qdk(n, u=2.0)
    classic = _classic(**CONVERGED, truncation_tol=1e-14)
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
    energy, _ = calculator.run(_dense_hamiltonian(h1, eri, 0.713), na, nb)
    from pyscf import fci

    reference = fci.direct_spin1.FCI().kernel(h1, eri, n, (na, nb))[0] + 0.713
    assert energy == pytest.approx(reference, abs=1e-9)
    _assert_rdms_match_fci(calculator, h1, eri, na, nb)


def test_classic_requires_successful_run_for_rdms():
    with pytest.raises(RuntimeError, match="has not produced"):
        _classic().make_rdms()


def test_classic_rdm_container_roundtrip(hubbard_qdk, tmp_path):
    from qdk_chemistry.data import Wavefunction

    calculator = _classic(sweeps=4)
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
