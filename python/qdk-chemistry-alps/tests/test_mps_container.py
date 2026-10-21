# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx tests/test_mps_container.py at 1adba871e:
# independent CAR permutation reference; classic final-state and QR gauge checks.
"""The optional MPS bridge preserves the final state, signs, RDMs and gauge."""
from itertools import product
from types import SimpleNamespace

import numpy as np
import pytest

from qdk_chemistry_alps.classic_dmrg_calculator import AlpsClassicDmrgCalculator


def blocked(tensors, alpha_first):
    n = len(tensors)
    result = np.zeros(4**n, dtype=complex)
    for physical in product(range(4), repeat=n):
        amplitude, created = np.ones((1, 1), dtype=complex), []
        for site, (tensor, local) in enumerate(zip(tensors, physical, strict=True)):
            amplitude = amplitude @ tensor[:, local, :]
            for spin in ((0, 1) if alpha_first else (1, 0)):
                if local & (1 << spin):
                    created.append(spin*n + site)
        inversions = sum(a > b for i, a in enumerate(created) for b in created[i+1:])
        result[sum(1 << mode for mode in created)] = (-1)**inversions*amplitude[0, 0]
    return result


def calculator(mode="mps", bond=16):
    c = AlpsClassicDmrgCalculator()
    for k, v in dict(wavefunction_container=mode, max_bond_dim=bond, sweeps=4,
                     ci_memory_limit_gb=0.001).items():
        c.settings().set(k, v)
    return c


@pytest.mark.classic_dmrg
def test_missing_capability_fails_before_solving(hubbard_qdk, monkeypatch):
    import qdk_chemistry.data as data
    from pyalps import classic_dmrg
    monkeypatch.delattr(data, "MPSContainer", raising=False)
    def forbidden(*args, **kwargs):
        raise AssertionError("must reject missing capability before native solve")
    monkeypatch.setattr(classic_dmrg, "solve_molecular", forbidden)
    with pytest.raises(RuntimeError, match="MPSContainer"):
        calculator().run(hubbard_qdk(2), 1, 1)


@pytest.mark.mps_qdk
@pytest.mark.classic_dmrg
@pytest.mark.parametrize(("n", "na", "nb", "bond"), [(2, 1, 1, 16), (3, 2, 1, 16), (4, 2, 2, 2)])
def test_same_final_state_and_serialization(hubbard_qdk, tmp_path, n, na, nb, bond):
    from qdk_chemistry.data import Wavefunction
    c = calculator(bond=bond)
    energy, w = c.run(hubbard_qdk(n), na, nb)
    assert len(c.last_result.ci_coefficients) == 0
    native = blocked(c.last_state.tensors, False)
    container = w.get_container()
    arrays = [s.to_dense().reshape(s.shape) for s in container.sites]
    actual = blocked(arrays, True)
    np.testing.assert_allclose(actual, native, atol=2e-12, rtol=0)
    assert container.orthogonality_center == 0
    for a in arrays[1:]:
        matrix = a.reshape(a.shape[0], -1)
        np.testing.assert_allclose(matrix @ matrix.conj().T, np.eye(a.shape[0]), atol=2e-12)
    dm1, dm2 = c.make_rdms()
    h1 = np.diag([-1.]*(n-1), 1)+np.diag([-1.]*(n-1), -1)
    erdm = np.einsum("pq,pq", h1, dm1)+2*np.einsum("iiii", dm2)
    assert abs(erdm-energy) < 1e-9
    aa, bb = w.get_active_one_rdm_spin_dependent()
    np.testing.assert_allclose([np.trace(aa), np.trace(bb)], [na, nb], atol=1e-9)
    np.testing.assert_allclose(w.get_active_two_rdm_spin_traced(), dm2.ravel(), atol=1e-12)
    path = tmp_path/'classic.wavefunction.h5'
    w.to_file(str(path), 'hdf5')
    loaded = Wavefunction.from_file(str(path), 'hdf5')
    assert loaded.content_hash() == w.content_hash()
    # Compare a second solve's CI vector after independently converting its creators.
    ci = calculator(mode="statevector", bond=bond)
    eci, wci = ci.run(hubbard_qdk(n), na, nb)
    reference = np.zeros(4**n, dtype=complex)
    for det, coefficient in zip(wci.get_active_determinants(), wci.get_coefficients(), strict=True):
        index = sum(int(bit) << q for q, bit in enumerate(det.to_bits(2*n)))
        reference[index] = coefficient
    assert abs(abs(np.vdot(reference, actual))-1) < 1e-9
    assert abs(energy-eci) < 1e-9


@pytest.mark.mps_qdk
@pytest.mark.parametrize("complex_values", [False, True])
def test_qr_and_basis_conversion_preserve_phase_scale_and_input(complex_values):
    from qdk_chemistry_alps.mps import _mps_sites
    rng = np.random.default_rng(174)
    arrays = [rng.normal(size=(1,4,3)), rng.normal(size=(3,4,1))]
    if complex_values:
        arrays = [a+1j*rng.normal(size=a.shape) for a in arrays]
    before = [a.copy() for a in arrays]
    state = SimpleNamespace(tensors=arrays)
    sites = _mps_sites(state)
    actual = blocked([s.to_dense().reshape(s.shape) for s in sites], True)
    np.testing.assert_allclose(actual, blocked(arrays, False), atol=1e-12, rtol=0)
    for a, b in zip(arrays, before, strict=True):
        np.testing.assert_array_equal(a,b)
