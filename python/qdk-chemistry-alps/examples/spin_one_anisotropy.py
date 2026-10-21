# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/spin_one_anisotropy.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 e073a448df4bb98694489930a9364f5dc74da64b98b66eba02bffe6445906cec.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Spin-1 anisotropy: changing local spin populations and bond correlations."""

import time

import numpy as np
from demo_support import export_operator
from pyalps import models
from demo_support import model_matrix
from pyalps.models import make_qubit_encoding
from pyalps import classic_dmrg as dmrg
from demo_support import finish, parser
from demo_support.spin_one import (
    COUPLINGS,
    FIELDS,
    model,
    qubit_matrix,
    reference_indices,
)

from qdk_chemistry_alps.interop import pauli_operator_to_qubit_operator


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout="constrained")
    D = data["anisotropy"]
    axes[0].plot(D, data["reference_energy"], color="#087e8b", label="Independent exact reference")
    axes[0].plot(D, data["dmrg_energy"], "+", color="#c44e52", label="ALPS DMRG")
    axes[0].plot(D, data["binary_energy"], "o", fillstyle="none", color="#d68c22", label="QDK binary physical block")
    axes[0].plot(D, data["gray_energy"], "x", color="#7956a8", label="QDK Gray physical block")
    axes[0].set(xlabel="Single-ion anisotropy D / J₀", ylabel="Ground-state energy / J₀", title="One model, two qubit encodings")
    axes[0].legend(fontsize=8)
    axes[1].plot(D, data["mean_sz_squared"], "o-", color="#087e8b", label="Mean ⟨(Sᶻ)²⟩")
    axes[1].plot(D, data["bond_sz_correlation"], "s-", color="#7956a8", label="⟨Sᶻ₀ Sᶻ₁⟩")
    axes[1].set(xlabel="Single-ion anisotropy D / J₀", ylabel="Expectation", title="Positive D favors the m = 0 state")
    axes[1].legend()
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    anisotropies = np.array([-1., 0., 1., 2., 4.])
    data = {"anisotropy": anisotropies, "reference_energy": [], "dmrg_energy": [],
            "binary_energy": [], "gray_energy": [], "mean_sz_squared": [], "bond_sz_correlation": [],
            "encoding_diagnostics": {}, "maximum_operator_error": 0.}
    for D in anisotropies:
        h, reference, axes = model(D)
        np.testing.assert_allclose(model_matrix(h), reference, atol=1e-13, rtol=0)
        energies, vectors = np.linalg.eigh(reference)
        psi = vectors[:, 0]
        spin = h.site_for(0)
        square = sum(z@z for z in axes[2])/4
        bond = axes[2][0] @ axes[2][1]
        total_sz = sum(axes[2])
        np.testing.assert_allclose(total_sz @ psi, 0., atol=1e-10)
        expression = sum((models.located(spin, "Sz", i)*models.located(spin, "Sz", i) for i in range(4)), start=models.TermSum())/4
        result = dmrg.solve(h, sector={"Sz": 0}, sweeps=dmrg.fixed_sweeps(6, 32),
                            expectation={"square": expression,
                                         "bond": models.located(spin, "Sz", 0)*models.located(spin, "Sz", 1)})
        np.testing.assert_allclose(result.energies[0], energies[0], atol=1e-9, rtol=0)
        expected_square = np.vdot(psi, square @ psi).real
        expected_bond = np.vdot(psi, bond @ psi).real
        np.testing.assert_allclose(result.measurements["square"].values[0,0], expected_square, atol=1e-8)
        np.testing.assert_allclose(result.measurements["bond"].values[0,0], expected_bond, atol=1e-8)
        data["reference_energy"].append(energies[0])
        data["dmrg_energy"].append(result.energies[0])
        data["mean_sz_squared"].append(expected_square)
        data["bond_sz_correlation"].append(expected_bond)
        for scheme in ("binary", "gray"):
            encoding = make_qubit_encoding(h, scheme=scheme)
            encoded = export_operator(h, encoding=encoding, tolerance=0.)
            qdk_h = pauli_operator_to_qubit_operator(encoded)
            matrix = qubit_matrix(qdk_h)
            indices = reference_indices(scheme)
            V = np.eye(256)[:, indices]
            error = np.max(abs(matrix @ V - V @ reference))
            assert error < 1e-12
            data["maximum_operator_error"] = max(data["maximum_operator_error"], error)
            mapped = encoding.encode_statevector(psi)
            np.testing.assert_allclose(mapped, V @ psi, atol=1e-14)
            assert encoding.physical_weight(mapped) > 1-1e-13
            observable = export_operator(h, observable=expression, encoding=encoding, tolerance=0.)
            mapped_square = qubit_matrix(pauli_operator_to_qubit_operator(observable))
            np.testing.assert_allclose(np.vdot(mapped, mapped_square @ mapped), expected_square, atol=1e-12)
            e = np.linalg.eigvalsh(matrix[np.ix_(indices, indices)])
            np.testing.assert_allclose(e, energies, atol=1e-12)
            # Trace completion is a unital positive map: the register's lowest
            # energy agrees too, but its ground eigenspace need not be physical.
            np.testing.assert_allclose(np.linalg.eigvalsh(matrix)[0], energies[0], atol=1e-12)
            data[scheme+"_energy"].append(e[0])
            data["encoding_diagnostics"][scheme] = {"layout": {"scheme": scheme, "site_dimensions": list(encoding.site_dimensions), "site_offsets": list(encoding.site_offsets)},
                "qubits": encoded.num_qubits, "pauli_terms": len(encoded.terms),
                "local_decompositions": encoded.local_decompositions, "local_cache_hits": encoded.local_cache_hits,
                "discarded_l1": encoded.discarded_l1}
    assert data["mean_sz_squared"][-1] < data["mean_sz_squared"][0] - .3
    print("Four spin-1 sites → eight qubits; binary and Gray physical spectra and observables agree.")
    print(f"Mean ⟨(Sz)²⟩: {data['mean_sz_squared'][0]:.4f} at D=-1 → {data['mean_sz_squared'][-1]:.4f} at D=4")
    print(f"Maximum mapped-operator error: {data['maximum_operator_error']:.3g}")
    print("QDK operator interoperability with ALPS DMRG; no quantum circuit is executed in this example.")
    data["model"] = {"sites": 4, "spin": 1, "boundary": "open", "J": COUPLINGS,
                     "K": .18, "h": FIELDS, "constant": .37, "units": "hbar=1, J0=1"}
    data["dmrg_controls"] = {"sector": {"Sz": 0}, "max_bond_dim": 32, "sweeps": 6, "lanczos_tolerance": 1e-10}
    finish(args, data, started, plot)


if __name__ == "__main__":
    main()
