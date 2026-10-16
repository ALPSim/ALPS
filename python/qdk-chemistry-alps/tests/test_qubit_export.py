# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_qubit_export.py. Changes: pyalps names;
# dense matrices come from the fock_matrix fixture.
"""Fermionic hardware export: ALPS Hamiltonian -> QubitOperator.

The strong check is spectral: the exported qubit operator must have the
same spectrum on the full 2^(2K) qubit space as the fermionic Hamiltonian
has on Fock space, constant term included. The second check pins our
conventions against qdk-chemistry's own mapper on identical integrals.
"""

from __future__ import annotations

import numpy as np
import pytest

from qdk_chemistry_alps import interop
from qdk_chemistry_alps.abinitio import build_abinitio_hamiltonian

ENCODINGS = ["jordan-wigner", "bravyi-kitaev", "bravyi-kitaev-tree", "parity"]


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


@pytest.fixture
def qubit_spectrum(fock_matrix):
    def spectrum(qubit_op) -> np.ndarray:
        H = interop.qubit_operator_to_hamiltonian(qubit_op)
        return np.sort(np.linalg.eigvalsh(fock_matrix(H)))

    return spectrum


@pytest.mark.parametrize("encoding", ENCODINGS)
def test_exported_spectrum_equals_the_fock_spectrum(
    encoding, hubbard_alps, fock_matrix, qubit_spectrum
):
    H = hubbard_alps(3, u=4.0, mu=0.7)
    fock = np.sort(np.linalg.eigvalsh(fock_matrix(H)))

    qo = interop.fermionic_hamiltonian_to_qubit_operator(H, encoding)
    assert qo.num_qubits == 6
    assert np.max(np.abs(qubit_spectrum(qo) - fock)) < 1e-12


def test_export_carries_encoding_metadata(hubbard_alps):
    qo = interop.fermionic_hamiltonian_to_qubit_operator(hubbard_alps(2))
    meta = qo.to_dict()
    assert meta["encoding"] == "jordan-wigner"
    assert meta["fermion_mode_order"] == "blocked"


def test_prebuilt_mapping_is_accepted(hubbard_alps, fock_matrix, qubit_spectrum):
    from qdk_chemistry.data import MajoranaMapping

    H = hubbard_alps(3)
    fock = np.sort(np.linalg.eigvalsh(fock_matrix(H)))
    mapping = MajoranaMapping.bravyi_kitaev(num_modes=6)

    qo = interop.fermionic_hamiltonian_to_qubit_operator(H, mapping)
    assert qo.to_dict()["encoding"] == "bravyi-kitaev"
    assert np.max(np.abs(qubit_spectrum(qo) - fock)) < 1e-12


def test_matches_qdk_own_mapper_on_identical_integrals(hubbard_qdk):
    """Our extraction round trip reproduces QdkQubitMapper term by term."""
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import MajoranaMapping

    n = 4
    qdk_h = hubbard_qdk(n)
    mapping = MajoranaMapping.jordan_wigner(num_modes=2 * n)

    reference = _terms(create("qubit_mapper").run(qdk_h, mapping))
    ecore = float(qdk_h.get_core_energy())
    if ecore:  # the mapper's Hamiltonian overload excludes the constant
        identity = "I" * (2 * n)
        reference[identity] = reference.get(identity, 0.0) + ecore

    h1 = np.asarray(qdk_h.get_one_body_integrals()[0])
    two = qdk_h.get_two_body_integrals()
    eri = np.asarray(two[0] if isinstance(two, (tuple, list)) else two)
    ours = _terms(
        interop.fermionic_hamiltonian_to_qubit_operator(
            build_abinitio_hamiltonian(h1, eri.reshape(n, n, n, n)), mapping
        )
    )

    assert set(ours) == set(reference)
    assert max(abs(ours[k] - reference[k]) for k in reference) < 1e-12


# --- error paths -----------------------------------------------------------

def test_unknown_encoding_lists_the_valid_ones(hubbard_alps):
    with pytest.raises(ValueError, match="unknown encoding"):
        interop.fermionic_hamiltonian_to_qubit_operator(hubbard_alps(2), "qubits-r-us")


def test_mapping_with_wrong_mode_count_is_rejected(hubbard_alps):
    from qdk_chemistry.data import MajoranaMapping

    with pytest.raises(ValueError, match="modes but the"):
        interop.fermionic_hamiltonian_to_qubit_operator(
            hubbard_alps(2), MajoranaMapping.jordan_wigner(num_modes=6)
        )


def test_tapered_mapping_points_at_the_qdk_mapper(hubbard_alps):
    from qdk_chemistry.data import MajoranaMapping, Symmetries

    mapping = MajoranaMapping.symmetry_conserving_bravyi_kitaev(
        num_modes=4, symmetries=Symmetries(1, 1)
    )
    if not mapping.tapering:
        pytest.skip("this qdk-chemistry build reports no tapering metadata")
    with pytest.raises(ValueError, match="tapering"):
        interop.fermionic_hamiltonian_to_qubit_operator(hubbard_alps(2), mapping)


def test_spin_hamiltonian_is_rejected():
    from pyalps import models

    with pytest.raises(ValueError, match="spinful-fermion sites"):
        interop.fermionic_hamiltonian_to_qubit_operator(
            models.heisenberg(models.lattice.chain(2, bc="open"), J=1.0)
        )


# --- placing a state into the qubit register -------------------------------
# Spectrum comparisons are invariant under any permutation of the basis, so
# they cannot see an endianness or parity error. These can.

def _gradient_hubbard(n: int, u: float = 4.0, slope: float = 0.6):
    """Hubbard chain with a site-energy gradient, so the density is not
    mirror-symmetric and a bit-reversed register is detectable."""
    h1 = np.zeros((n, n))
    for i in range(n):
        h1[i, i] = slope * (i - (n - 1) / 2)
    for i in range(n - 1):
        h1[i, i + 1] = h1[i + 1, i] = -1.0
    eri = np.zeros((n, n, n, n))
    for i in range(n):
        eri[i, i, i, i] = u
    return h1, eri


def _mapped_state(n: int, n_alpha: int, n_beta: int):
    from qdk_chemistry_alps.abinitio import solve_ground_state

    h1, eri = _gradient_hubbard(n)
    H = build_abinitio_hamiltonian(h1, eri)
    result = solve_ground_state(H, n_alpha, n_beta)
    psi = interop.ci_vector_to_qubit_statevector(
        result.problem.states, np.asarray(result.eigenvectors)[:, 0], n
    )
    return H, result, psi


@pytest.fixture
def qubit_matrix(fock_matrix):
    def matrix(qubit_op, n_qubits: int) -> np.ndarray:
        """Dense operator, reindexed from ALPS's ordering to little-endian."""
        dense = fock_matrix(interop.qubit_operator_to_hamiltonian(qubit_op))
        order = np.array(
            [int(format(i, f"0{n_qubits}b")[::-1], 2) for i in range(2**n_qubits)]
        )
        return dense[np.ix_(order, order)]

    return matrix


def test_mapped_state_reproduces_the_alps_energy(qubit_matrix):
    """The decisive convention check: same state, same operator, same energy."""
    n = 4
    H, result, psi = _mapped_state(n, 2, 2)
    assert abs(np.linalg.norm(psi) - 1.0) < 1e-10

    qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H)
    matrix = qubit_matrix(qubit_op, 2 * n)
    energy = float(np.real(np.vdot(psi, matrix @ psi)))
    assert abs(energy - float(result.eigenvalues[0])) < 1e-9


def test_a_bit_reversed_register_is_detected(qubit_matrix):
    """Negative control: the check above must fail on a mirrored mapping.

    Without this, a passing test proves nothing -- a spectrum comparison
    passes under bit reversal too.
    """
    n = 4
    H, result, psi = _mapped_state(n, 2, 2)
    qubit_op = interop.fermionic_hamiltonian_to_qubit_operator(H)
    matrix = qubit_matrix(qubit_op, 2 * n)

    reversed_index = np.array(
        [int(format(i, f"0{2 * n}b")[::-1], 2) for i in range(2 ** (2 * n))]
    )
    mirrored = psi[reversed_index]
    energy = float(np.real(np.vdot(mirrored, matrix @ mirrored)))
    assert abs(energy - float(result.eigenvalues[0])) > 1e-3, (
        "a bit-reversed state gave the right energy, so this check is blind "
        "to endianness -- the model is probably mirror-symmetric"
    )


def test_mapped_state_has_definite_particle_number():
    """Every populated amplitude carries the electron count of the sector.

    This is what makes particle-number post-selection valid on a device.
    """
    n, n_alpha, n_beta = 4, 2, 2
    _, _, psi = _mapped_state(n, n_alpha, n_beta)
    populated = np.nonzero(np.abs(psi) > 1e-12)[0]
    assert populated.size
    for index in populated:
        assert bin(int(index)).count("1") == n_alpha + n_beta


def test_state_mapping_validates_its_inputs():
    states = np.zeros((3, 4), dtype=np.int64)
    with pytest.raises(ValueError, match="basis states but"):
        interop.ci_vector_to_qubit_statevector(states, np.zeros(2), 4)
    with pytest.raises(ValueError, match="span"):
        interop.ci_vector_to_qubit_statevector(states, np.zeros(3), 5)
