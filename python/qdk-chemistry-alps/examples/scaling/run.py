# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/scaling/run.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 2e18f89ce3f7c39b0dd62fb4380fbf297dce7fb3fe77817259512cb5345d6e51.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Measured ALPS/QDK scaling examples; see README.md for the meeting runs.

Large calculations require QDK MPSContainer or RdmContainer. No CI vector is
reconstructed. CASCI/FCI references are limited to <= 50,000 determinants.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import subprocess
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
from qdk_chemistry.algorithms import create
from qdk_chemistry.data import (
    CanonicalFourCenterHamiltonianContainer,
    Hamiltonian,
    LatticeGraph,
    ModelOrbitals,
)
from qdk_chemistry.utils.model_hamiltonians import create_hubbard_hamiltonian

import qdk_chemistry_alps  # noqa: F401 -- register the calculator


def problem(args):
    n = args.sites
    if args.case == "hubbard":
        ham = create_hubbard_hamiltonian(
            LatticeGraph.chain(n, periodic=False), epsilon=0.0, t=1.0, U=args.u
        )
        return ham, {"name": f"Hubbard {n}", "boundary": "open", "t": 1.0,
                     "U": args.u, "energy_unit": "t"}

    from pyscf import ao2mo, gto, scf
    from scipy.linalg import fractional_matrix_power

    # Full STO-3G space: every electron and orbital is active. Symmetric
    # Lowdin orthogonalization keeps orbitals ordered along the chain.
    mol = gto.M(atom=[("H", (0, 0, i * args.spacing)) for i in range(n)],
                basis="sto-3g", unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-11
    mf.max_cycle = 200
    mf.kernel()
    used_newton = not mf.converged
    if used_newton:
        # Long, stretched chains can stall DIIS near orbital degeneracies.
        # Newton SCF supplies the RHF reference; the active orbitals below
        # still come from the overlap matrix, independently of this reference.
        mf = mf.newton()
        mf.max_cycle = 100
        mf.kernel()
    if not mf.converged:
        raise RuntimeError("RHF did not converge")
    c = fractional_matrix_power(mf.get_ovlp(), -0.5).real
    np.testing.assert_allclose(c.T @ mf.get_ovlp() @ c, np.eye(n), atol=1e-10)
    h1 = c.T @ mf.get_hcore() @ c
    eri = ao2mo.restore(1, ao2mo.kernel(mol, c), n)
    # Restore exact real-integral permutation symmetries after floating-point
    # AO transformations. Otherwise roundoff sends the generic Hermiticity
    # validator through its bounded, much more expensive word expansion.
    symmetric_h1 = 0.5 * (h1 + h1.T)
    symmetric_eri = ao2mo.restore(1, ao2mo.restore(8, eri, n), n)
    symmetry_roundoff = max(float(np.max(abs(h1 - symmetric_h1))),
                            float(np.max(abs(eri - symmetric_eri))))
    if symmetry_roundoff > 1e-12:
        raise ValueError("integrals violate real permutation symmetry beyond roundoff")
    h1, eri = symmetric_h1, symmetric_eri
    ham = Hamiltonian(CanonicalFourCenterHamiltonianContainer(
        h1, eri.reshape(-1), ModelOrbitals(n), mol.energy_nuc(), h1
    ))
    return ham, {"name": f"H{n} / STO-3G", "spacing_angstrom": args.spacing,
                 "basis": "STO-3G", "orbital_order": "symmetric Lowdin AOs, chain order",
                 "rhf_energy": float(mf.e_tot), "nuclear_energy": float(mol.energy_nuc()),
                 "rhf_newton_fallback": used_newton,
                 "integral_symmetry_roundoff_removed": symmetry_roundoff,
                 "energy_unit": "Eh", "geometry_angstrom": mol.atom_coords(unit="Angstrom").tolist()}


def integrals(ham, n):
    h1 = np.asarray(ham.get_one_body_integrals()[0])
    two = ham.get_two_body_integrals()
    eri = np.asarray(two[0] if isinstance(two, (tuple, list)) else two).reshape((n,) * 4)
    return h1, eri, float(ham.get_core_energy())


def peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss * 1024)


def provenance():
    import pyalps
    import qdk_chemistry

    import qdk_chemistry_alps.classic_dmrg_calculator as adapter

    repo = Path(__file__).resolve().parents[4]
    revision = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=False).stdout.strip()
    versions = {p: importlib.metadata.version(p)
                for p in ("pyalps", "qdk-chemistry", "numpy", "scipy", "pyscf")}
    return {"generated_utc": datetime.datetime.now(datetime.UTC).isoformat(),
            "python": sys.version, "executable": sys.executable, "platform": platform.platform(),
            "versions": versions, "alps_module": pyalps.__file__, "qdk_module": qdk_chemistry.__file__,
            "adapter_module": adapter.__file__, "alps_repository_revision": revision,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "adapter_sha256": hashlib.sha256(Path(adapter.__file__).read_bytes()).hexdigest(),
            "thread_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
            "memory_note": "process cumulative peak RSS, not per-solver memory; includes imports and integrals"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("hubbard", "hydrogen"), required=True)
    parser.add_argument("--sites", type=int, default=8)
    parser.add_argument("--bonds", type=int, nargs="+", default=[32, 64])
    parser.add_argument("--sweeps", type=int, default=12)
    parser.add_argument("--container", choices=("mps", "rdm"), default="mps")
    parser.add_argument("--u", type=float, default=8.0)
    parser.add_argument("--spacing", type=float, default=1.8)
    parser.add_argument("--integral-cutoff", type=float, default=1e-12)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.sites < 2 or args.sites % 2 or min(args.bonds) < 1 or args.sweeps < 2:
        parser.error("use an even site count >= 2, positive bonds and >= 2 sweeps")
    # Fail before spending time on a solve with an unsupported QDK release.
    if args.container == "mps":
        from qdk_chemistry_alps.mps import require_mps_container
        require_mps_container()
    else:
        from qdk_chemistry.data import RdmContainer  # noqa: F401

    started = time.perf_counter()
    ham, description = problem(args)
    n = args.sites
    h1, eri, core = integrals(ham, n)
    if not np.isfinite(args.integral_cutoff) or args.integral_cutoff < 0:
        parser.error("integral cutoff must be finite and nonnegative")
    full_h1, full_eri = h1, eri
    h1 = np.where(abs(h1) >= args.integral_cutoff, h1, 0.0)
    eri = np.where(abs(eri) >= args.integral_cutoff, eri, 0.0)
    # Each fermionic monomial has operator norm <= 1. Two spin channels
    # in h1 and four channels times 1/2 in eri give this conservative bound
    # on both any state's energy shift and the exact ground-energy shift.
    screening_bound = 2 * float(abs(full_h1 - h1).sum() + abs(full_eri - eri).sum())
    ham = Hamiltonian(CanonicalFourCenterHamiltonianContainer(
        h1, eri.reshape(-1), ham.get_orbitals(), core, h1
    ))
    determinants = comb(n, n // 2) ** 2
    data = {"problem": description | {"orbitals": n, "electrons": n, "n_alpha": n // 2,
                                      "n_beta": n // 2, "determinants": determinants,
                                      "ci_float64_vector_bytes": 8 * determinants},
            "controls": {"sweeps": args.sweeps,  "eigensolver_tol": 1e-10,
                         "truncation_tol": 0.0, "integral_cutoff": args.integral_cutoff},
            "screening_energy_bound": screening_bound,
            "preparation_seconds": time.perf_counter() - started,
            "provenance": provenance(), "dmrg": []}

    def save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(data, indent=2) + "\n")

    print(f"{description['name']}: CAS({n},{n}), {determinants:,} determinants; "
          f"one real CI vector = {8 * determinants / 1e12:.6g} TB", flush=True)
    if determinants <= 50_000:
        calc = create("multi_configuration_calculator", "macis_cas")
        calc.settings().set("calculate_one_rdm", True)
        calc.settings().set("calculate_two_rdm", True)
        calc.settings().set("ci_residual_tolerance", 1e-10)
        t0 = time.perf_counter()
        exact, _ = calc.run(ham, n // 2, n // 2)
        data["casci"] = {"status": "measured", "energy": exact,
                         "seconds_including_rdms": time.perf_counter() - t0}
        from pyscf import fci
        solver = fci.direct_spin1.FCI()
        solver.conv_tol = 1e-12
        e_fci, ci = solver.kernel(h1, eri, n, (n // 2, n // 2), ecore=core)
        ref_dm1, _ = solver.make_rdm12(ci, n, (n // 2, n // 2))
        data["pyscf_fci"] = {"energy": e_fci}
        np.testing.assert_allclose(exact, e_fci, atol=1e-8, rtol=0)
        print(f"MACIS CASCI: {exact:.12f}; PySCF FCI: {e_fci:.12f}", flush=True)
    else:
        data["casci"] = {"status": "not_run", "reason": "determinant count exceeds the 50,000 reference guard; vector memory is an analytical lower bound, not a measured failure"}
    save()

    for bond in args.bonds:
        calc = create("multi_configuration_calculator", "alps_classic_dmrg")
        for key, value in {"wavefunction_container": args.container, "max_bond_dim": bond,
                           "sweeps": args.sweeps, "truncation_tol": 0.0,
                           "integral_cutoff": args.integral_cutoff}.items():
            calc.settings().set(key, value)
        t0 = time.perf_counter()
        energy, wfn = calc.run(ham, n // 2, n // 2)
        elapsed = time.perf_counter() - t0
        container = wfn.get_container()
        dm1 = np.asarray(wfn.get_active_one_rdm_spin_traced())
        dm2 = np.asarray(wfn.get_active_two_rdm_spin_traced()).reshape((n,) * 4)
        reconstructed = float(np.einsum("pq,pq", h1, dm1)
                              + 0.5 * np.einsum("pqrs,pqrs", eri, dm2) + core)
        full_energy = float(np.einsum("pq,pq", full_h1, dm1)
                            + 0.5 * np.einsum("pqrs,pqrs", full_eri, dm2) + core)
        occupations = np.linalg.eigvalsh(dm1)[::-1]
        np.testing.assert_allclose(reconstructed, energy, atol=2e-7, rtol=0)
        assert abs(full_energy - energy) <= screening_bound + 2e-7
        np.testing.assert_allclose(np.trace(dm1), n, atol=1e-8, rtol=0)
        assert occupations.min() >= -1e-8 and occupations.max() <= 2 + 1e-8
        assert wfn.get_container_type() == args.container
        history = np.asarray([it.energy for it in calc.last_result.iterations])
        tensors = calc.last_state.tensors if calc.last_state is not None else []
        record = {"bond_requested": bond, "bond_actual": max(calc.last_state.bond_dimensions) if calc.last_state is not None else None,
                  "energy": energy, "seconds_including_rdms": elapsed,
                  "energy_full_integrals_from_rdms": full_energy,
                  "process_peak_rss_bytes": peak_rss_bytes(),
                  "dense_mps_export_bytes": sum(t.nbytes for t in tensors),
                  "dense_mps_parameters": sum(t.size for t in tensors),
                  "returned_container": wfn.get_container_type(),
                  "rdm_energy_error": abs(reconstructed - energy),
                  "rdm_trace_error": float(abs(np.trace(dm1) - n)),
                  "natural_occupations": occupations.tolist(),
                  "double_occupancies": (np.einsum("iiii->i", dm2) / 2).tolist(),
                  "iteration_energies_electronic": history.tolist(),
                  "last_iteration_energy_change": float(abs(history[-1] - history[-2])),
                  "maximum_iteration_discarded_weight": max(it.truncation_error for it in calc.last_result.iterations),
                  "endpoint_discarded_weight": calc.last_result.final_truncation_error,
                  "convergence_note": "adjacent iteration energies are not a completed-sweep convergence criterion"}
        if "pyscf_fci" in data:
            record["error_vs_casci"] = energy - data["casci"]["energy"]
            record["one_rdm_max_error_vs_fci"] = float(np.max(np.abs(dm1 - ref_dm1)))
            assert energy >= data["casci"]["energy"] - 1e-8
        data["dmrg"].append(record)
        save()
        print(f"D={bond}: E={energy:.12f}, {elapsed:.2f} s incl. RDMs; "
              f"last iteration dE={record['last_iteration_energy_change']:.2e}; "
              f"RDM energy error={record['rdm_energy_error']:.2e}", flush=True)
        del tensors, calc, wfn, container, dm1, dm2
    print(f"Saved {args.output}", flush=True)


if __name__ == "__main__":
    main()
