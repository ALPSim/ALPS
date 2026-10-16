# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_pauli_interop.py. Changes: pyalps
# names (models.export_pauli, compile_operator and physical_basis for
# observables, pyalps.ed for the physical matrices).
"""Container-only conversion of the native general Pauli export."""

import numpy as np
import pytest
from pyalps import ed, models

from qdk_chemistry_alps import interop


def _qdk_register_matrix(qdk, fock_matrix, n_qubits):
    """Dense QubitOperator matrix in QDK register order (qubit 0 = bit 0).

    The spin importer builds product order with site 0 most significant;
    undo that bit reversal.
    """
    reverse = [int(f"{i:0{n_qubits}b}"[::-1], 2) for i in range(2**n_qubits)]
    dense = fock_matrix(interop.qubit_operator_to_hamiltonian(qdk))
    return dense[np.ix_(reverse, reverse)]


def test_pairing_export_uses_native_layout_without_integrals(monkeypatch, fock_matrix):
    def forbidden(*args, **kwargs):
        raise AssertionError("general export must not allocate chemistry integrals")

    monkeypatch.setattr(interop, "hamiltonian_to_integrals", forbidden)
    site = models.spinless_fermion()

    def op(name, i):
        return models.located(site, name, i)

    h = models.Hamiltonian(models.lattice.chain(3, bc="open"), site)
    h.add((.2 + .13j) * (op("c", 2) * op("c", 0)), add_hc=True)
    h.add(.31 * op("n", 1))
    exported = models.export_pauli(h)
    qdk = interop.pauli_operator_to_qubit_operator(exported)
    data = qdk.to_dict()
    actual = dict(zip(data["pauli_strings"], np.asarray(data["coefficients"]["real"])
                      + 1j * np.asarray(data["coefficients"]["imag"]), strict=True))
    assert actual == dict(exported.terms)
    # The explicit native state map places each physical product state.
    problem = ed.prepare_product(h)
    indices = [exported.state_index(list(row)) for row in np.asarray(problem.states)]
    register = _qdk_register_matrix(qdk, fock_matrix, 3)
    np.testing.assert_allclose(register[np.ix_(indices, indices)],
                               problem.hamiltonian.to_scipy().toarray(), atol=1e-14)
    assert data.get("fermion_mode_order") is None


def test_empty_and_nonhermitian_operator_conversion():
    site = models.spinless_fermion()
    h = models.Hamiltonian(models.lattice.chain(2, bc="open"), site)
    zero = interop.pauli_operator_to_qubit_operator(models.export_pauli(h))
    assert zero.num_qubits == 2
    assert zero.to_dict()["pauli_strings"] == ["II"]
    odd = models.export_pauli(
        models.compile_operator(models.located(site, "c", 1), [site, site]),
        models.physical_basis(h),
        require_hermitian=False, require_even=False, require_conservation=False,
    )
    data = interop.pauli_operator_to_qubit_operator(odd).to_dict()
    assert np.max(np.abs(data["coefficients"]["imag"])) == .5


@pytest.mark.parametrize("scheme", ["binary", "gray"])
def test_spin_one_state_and_observable_use_the_same_native_encoding(scheme, fock_matrix):
    spin = models.spin(1.0)

    def op(name, i):
        return models.located(spin, name, i)

    h = models.Hamiltonian(models.lattice.chain(2, bc="open"), spin)
    h.add(.7 * (op("Sy", 0) * op("Sx", 1)) + .23 * (op("Sz", 1) * op("Sz", 1))
          + 1.8 * op("I", 0))
    encoding = models.make_qubit_encoding(h, scheme=scheme)
    y = np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]) / np.sqrt(2)
    x = abs(y)
    expected = .7 * np.kron(y, x) + .23 * np.kron(np.eye(3), np.diag([1, 0, 1])) + 1.8 * np.eye(9)
    # Build the isometry independently: site 0 in the low register bits.
    codes = [0, 1, 2] if scheme == "binary" else [0, 1, 3]
    indices = [codes[a] + 4 * codes[b] for a in range(3) for b in range(3)]
    V = np.eye(16)[:, indices]
    observable = models.export_pauli(
        models.compile_operator(op("Sp", 1), [spin, spin]), models.physical_basis(h),
        encoding=encoding, tolerance=0.0,
        require_hermitian=False, require_conservation=False,
    )
    for exported, physical in [
        (models.export_pauli(h, encoding=encoding, tolerance=0.0), expected),
        (observable, np.kron(np.eye(3), x + 1j * y)),
    ]:
        qdk = interop.pauli_operator_to_qubit_operator(exported)
        matrix = _qdk_register_matrix(qdk, fock_matrix, 4)
        np.testing.assert_allclose(matrix @ V, V @ physical, atol=1e-13)
        assert qdk.num_qubits == 4
    state = np.arange(9) + 1j * np.arange(9)[::-1]
    state /= np.linalg.norm(state)
    np.testing.assert_array_equal(encoding.encode_statevector(state), V @ state)
