# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, dmrg_calculator.py at 1adba871e,
# Source worktree snapshot SHA-256 841ebae37ea47f21b250c354ff2fbdfbac8f36daab3fd394121d267d3de62e31.
# _mps_sites and _mps_particle_counts. Uses classic MatrixProductState tensors
# and lossless QR sweeps to supply the center-zero gauge required by PR #540.
"""Optional QDK MPS export without a determinant expansion."""
import numpy as np


def require_mps_container():
    try:
        from qdk_chemistry.data import MPSContainer, MPSSite  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("wavefunction_container='mps' requires QDK MPSContainer and MPSSite "
                           "(PR #751); use 'statevector' for small problems or 'rdm' "
                           "with a QDK build providing RdmContainer") from exc
    return MPSContainer


def _right_canonical_tensors(state):
    """Move the center to zero with exact QR factorizations, never truncate."""
    tensors = [np.array(a, copy=True) for a in state.tensors]
    for i in range(len(tensors)-1, 0, -1):
        left, physical, right = tensors[i].shape
        q, r = np.linalg.qr(tensors[i].reshape(left, physical*right).conj().T, mode="reduced")
        tensors[i] = q.conj().T.reshape(q.shape[1], physical, right)
        tensors[i-1] = np.tensordot(tensors[i-1], r.conj().T, axes=(2, 0))
    return tensors


def _mps_sites(state):
    """Copy ALPS orbital tensors into QDK's alpha-before-beta local basis.

    Dense export loses auxiliary charge labels, so each site uses one trivial
    symmetry block. It does not expand the many-body wavefunction; QR can remove
    structurally redundant bond dimensions without truncating the state. ALPS owns the original state and its symmetry blocks.
    """
    from qdk_chemistry.data import MPSSite
    from qdk_chemistry.data import symmetry as sym

    product, label = sym.SymmetryProduct([]), sym.SymmetryLabel([])
    sites = []
    for values in _right_canonical_tensors(state):
        left, physical, right = values.shape
        if physical != 4:
            raise ValueError("the chemistry MPS bridge requires SpinfulFermion sites")
        # ALPS |2> = c†down c†up|0>; QDK |2> = c†up c†down|0>.
        # Only this local phase belongs here. Consumers handle the permutation
        # from site-interleaved modes to blocked alpha-then-beta modes.
        values[:, 3, :] *= -1
        tensor_type = (sym.SymmetryBlockedTensorRank3Complex
                       if np.iscomplexobj(values) else sym.SymmetryBlockedTensorRank3)
        tensor = tensor_type(
            [product] * 3,
            [{label: left}, {label: physical}, {label: right}],
            [((label,) * 3, values.reshape(left * physical, right))],
        )
        sites.append(MPSSite(tensor, [label], [label], [label]))
    return sites


def _mps_particle_counts(orbitals, n_alpha, n_beta):
    """Pack active and total counts using QDK's orbital symmetry convention."""
    from qdk_chemistry.data import symmetry as sym

    alpha = sym.SymmetryLabel([sym.axes.alpha()])
    beta = sym.SymmetryLabel([sym.axes.beta()])
    trivial = sym.SymmetryLabel([])
    spinful = orbitals.symmetries().has_axis(sym.AxisName.Spin)
    product = sym.SymmetryProduct([sym.axes.spin(1, False)] if spinful else [])

    def pack(a, b):
        blocks = [((alpha,), a), ((beta,), b)] if spinful else [((trivial,), a + b)]
        return sym.SymmetryBlockedScalarCount([product], blocks)

    inactive = orbitals.inactive_indices()
    totals = [n_alpha, n_beta]
    if inactive is not None:
        for i, label in enumerate((alpha, beta)):
            if not inactive.symmetries().has_axis(sym.AxisName.Spin):
                label = trivial
            if inactive.has(label):
                totals[i] += len(inactive.indices(label))
    return pack(n_alpha, n_beta), pack(*totals)



def mps_container(state, orbitals, n_alpha, n_beta, **rdms):
    """Copy the same state and its RDMs into Microsoft's optional container."""
    cls = require_mps_container()
    active, total = _mps_particle_counts(orbitals, n_alpha, n_beta)
    return cls(_mps_sites(state), orbitals, orthogonality_center=0,
               active_num_particles=active, total_num_particles=total, **rdms)
