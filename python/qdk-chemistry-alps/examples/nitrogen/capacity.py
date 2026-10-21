# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/nitrogen/capacity.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 7f0d7e0fd25f254b4ef0864581d0dc287b0f06ea2315c65fd8c07e2c01909641.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Frozen-core N2/cc-pVDZ: CAS(10e,26o), without a CI vector allocation."""
import argparse
import hashlib
import json
import tempfile
import time
from math import comb

import numpy as np

from run import rss, save



def solve_and_export(args, data):
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import CanonicalFourCenterHamiltonianContainer, Hamiltonian, ModelOrbitals, Wavefunction
    import qdk_chemistry_alps  # noqa: F401
    ham = Hamiltonian(CanonicalFourCenterHamiltonianContainer(
        data["h1"], data["eri"].ravel(), ModelOrbitals(26), float(data["core"]), data["h1"]))
    calc = create("multi_configuration_calculator", "alps_classic_dmrg",
                  wavefunction_container="mps", max_bond_dim=args.bond,
                  sweeps=args.sweeps, truncation_tol=0.0, convergence_tol=1e-10)
    start = time.perf_counter()
    energy, wfn = calc.run(ham, 5, 5)
    dm1, dm2 = calc.make_rdms()
    erdm = float(np.einsum("pq,pq", data["h1"], dm1)
                 + .5*np.einsum("pqrs,pqrs", data["eri"], dm2)+data["core"])
    np.testing.assert_allclose(erdm, energy, atol=1e-7, rtol=0)
    np.testing.assert_allclose(np.trace(dm1), 10, atol=1e-7, rtol=0)
    np.testing.assert_allclose(np.einsum("ppqq",dm2), 90, atol=1e-6, rtol=0)
    path=args.folder/f"d{args.bond}.wavefunction.h5"
    wfn.to_file(str(path), "hdf5")
    assert Wavefunction.from_file(str(path), "hdf5").content_hash()==wfn.content_hash()
    return {"energy": energy, "bond_requested": args.bond,
            "bond_actual": wfn.get_container().max_bond_dimension,
            "seconds_including_rdms_and_export": time.perf_counter()-start,
            "rdm_energy_error": abs(erdm-energy), "spin_squared": float(calc.spin_spin_correlations().sum()),
            "natural_occupations": np.linalg.eigvalsh(dm1)[::-1].tolist(),
            "mps_file_bytes": path.stat().st_size, "wavefunction_hash": wfn.content_hash(),
            "iterations": [{"index": it.index, "direction": str(it.direction), "energy": it.energy,
                            "discarded_weight": it.truncation_error} for it in calc.last_result.iterations],
            "note": "Independent fixed-state solve; no MPO compilation, MPS restart or one-site polishing API"}

def main():
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "block2", "alps"])
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--distance", type=float, default=2.0)
    parser.add_argument("--bond", type=int, default=64)
    parser.add_argument("--sweeps", type=int, default=24)
    args = parser.parse_args()
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.sweeps < 2 or args.bond < 1:
        parser.error("use at least two sweeps and a positive bond dimension")
    args.folder.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    path = args.folder / "integrals.npz"
    if args.stage == "prepare":
        from pyscf import ao2mo, gto, mcscf, scf
        mol = gto.M(atom=f"N 0 0 0; N 0 0 {args.distance}", basis="cc-pvdz", verbose=0,
                    max_memory=1500, unit="Angstrom", symmetry="D2h")
        mf = scf.RHF(mol).run(conv_tol=1e-12)
        assert mf.converged
        cas = mcscf.CASCI(mf, 26, 10)  # Integral transformation only; never kernel().
        h1, core = cas.get_h1eff()
        h1 = (h1+h1.T)/2
        eri = ao2mo.restore(1, ao2mo.restore(8, cas.get_h2eff(), 26), 26)
        # Exchange-weighted Fiedler ordering, the same interaction matrix as
        # Block2's orbital_reordering method. Apply explicitly to both solvers.
        interaction = abs(np.einsum("ijji->ij", eri)) + 1e-7*abs(h1)
        _, modes = np.linalg.eigh(np.diag(interaction.sum(axis=1))-interaction)
        order = np.argsort(modes[:, 1])
        h1, eri = h1[np.ix_(order, order)], eri[np.ix_(order, order, order, order)]
        cutoff = 1e-10
        hs, es = np.where(abs(h1) >= cutoff, h1, 0), np.where(abs(eri) >= cutoff, eri, 0)
        bound = 2*(abs(h1-hs).sum()+abs(eri-es).sum())
        np.savez_compressed(path, h1=hs, eri=es, core=core, mo_coeff=mf.mo_coeff[:, 2+order], order=order)
        record = {"distance_angstrom": args.distance, "basis": "cc-pVDZ", "electrons": 10,
                  "orbitals": 26, "frozen_core_electrons": 4,
                  "orbital_basis": "canonical RHF, D2h; exchange-weighted Fiedler ordering",
                  "orbital_permutation": order.tolist(),
                  "rhf_energy": mf.e_tot, "core_energy": float(core),
                  "determinants": comb(26, 5)**2, "one_float64_ci_vector_bytes": 8*comb(26, 5)**2,
                  "integral_cutoff": cutoff, "screening_energy_bound_Eh": float(bound),
                  "integrals_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "casci": "Not run: 4.327 billion determinants; the workflow avoids CI reconstruction"}
    else:
        data = np.load(path)
        if args.stage == "block2":
            from pyblock2.driver.core import DMRGDriver, SymmetryTypes
            with tempfile.TemporaryDirectory(prefix="n2-large-block2-") as scratch:
                driver = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SZ,
                                    n_threads=1, stack_mem=1_000_000_000)
                driver.initialize_system(n_sites=26, n_elec=10, spin=0)
                mpo = driver.get_qc_mpo(h1e=data["h1"], g2e=data["eri"], ecore=float(data["core"]), iprint=0)
                ket = driver.get_random_mps(tag="N2", bond_dim=args.bond, nroots=1)
                energy = driver.dmrg(mpo, ket, n_sweeps=args.sweeps, bond_dims=[args.bond],
                                     noises=[1e-5]*4+[1e-7]*4+[0], thrds=[1e-12], tol=1e-10, iprint=1)
                _, discarded, history = driver.get_dmrg_results()
                history = history.ravel()
                dm1 = np.asarray(driver.get_1pdm(ket, iprint=0)).sum(axis=0)
                record = {"energy": energy, "bond": args.bond, "sweeps_max": args.sweeps,
                          "sweep_energies": history.tolist(), "sweeps_completed": len(history),
                          "last_sweep_delta_Eh": float(abs(history[-1]-history[-2])),
                          "last_discarded_weight": float(discarded.ravel()[-1]),
                          "natural_occupations": np.linalg.eigvalsh(dm1)[::-1].tolist()}
        else:
            record = solve_and_export(args, data)
    record.update(seconds=time.perf_counter()-start, process_peak_rss_bytes=rss(),
                  script_sha256=script_hash,
                  integrals_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    suffix = "" if args.stage == "prepare" else f"-d{args.bond}"
    save(args.folder / f"{args.stage}{suffix}.json", record)
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
