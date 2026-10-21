# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/nitrogen/run.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 6ad6294b1c7b2f186df96abe15f7cff40d0f900aeed81235ed83e338b489c197.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""N2 bond breaking and the ALPS -> Microsoft MPS -> Q# workflow.

Run each stage in a fresh, memory-guarded process; see README.md.
Dense statevectors below are validation references for this 12-qubit example
only. Neither the ALPS MPS export nor the circuit constructor uses them.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def sector_indices(n=6, ne=3):
    from pyscf.fci import cistring
    strings = cistring.make_strings(range(n), ne)
    return (strings[:, None] | (strings[None, :] << n)).ravel()


def prepare(args):
    from pyscf import ao2mo, fci, gto, mcscf, scf
    from pyscf.mcscf import avas

    mol = gto.M(atom=f"N 0 0 0; N 0 0 {args.distance}", basis="cc-pvdz",
                unit="Angstrom", spin=0, charge=0, verbose=0, max_memory=1500)
    mf = scf.RHF(mol).run(conv_tol=1e-12)
    assert mf.converged
    # Exclude four doubly occupied orbitals (1s and 2s pairs) from the
    # initial AO projection. AVAS then selects the six 2p bonding/antibonding
    # orbitals. CASSCF optimizes inactive/active/virtual orbital rotations.
    n, ne, mo = avas.avas(mf, ["N 2p"], ncore=4, threshold=0.2)
    assert (n, ne) == (6, 6), (n, ne)
    mc = mcscf.CASSCF(mf, n, ne)
    mc.conv_tol, mc.conv_tol_grad, mc.max_cycle_macro = 1e-10, 1e-6, 100
    mc.fcisolver.conv_tol = 1e-12
    mc.kernel(mo)
    newton_fallback = not mc.converged
    if newton_fallback:
        mc = mc.newton()
        mc.max_cycle_macro = 100
        mc.kernel(mc.mo_coeff, mc.ci)
    assert mc.converged
    h1, core = mc.get_h1eff()
    h1 = (h1 + h1.T) / 2
    eri = ao2mo.restore(1, ao2mo.restore(8, mc.get_h2eff(), n), n)
    # Keep each bonding/antibonding pair adjacent. Match by the exchange
    # integrals, without using a CI vector or fitting to the DMRG result.
    from scipy.optimize import linear_sum_assignment
    _, partners = linear_sum_assignment(-abs(np.einsum("ijji->ij", eri)[:3, 3:]))
    order = np.array([v for i, j in enumerate(partners) for v in (i, 3+int(j))])
    h1, eri = h1[np.ix_(order, order)], eri[np.ix_(order, order, order, order)]
    mo_ordered = mc.mo_coeff.copy()
    mo_ordered[:, 4:10] = mc.mo_coeff[:, 4+order]
    solver = fci.direct_spin1.FCI()
    solver.conv_tol = 1e-13
    efci, ci = solver.kernel(h1, eri, n, (3, 3), ecore=core)
    dm1, dm2 = solver.make_rdm12(ci, n, (3, 3))
    np.testing.assert_allclose(efci, mc.e_tot, atol=1e-9, rtol=0)
    vector = np.zeros(4**n, dtype=complex)
    vector[sector_indices()] = ci.ravel()
    bonding_mask = sum(1 << i for i, old in enumerate(order) if old < 3)
    determinant_probability = abs(vector[bonding_mask | (bonding_mask << n)])**2
    # Actual overlap with the original RHF determinant, including the
    # optimized inactive orbitals. This is not just the largest CI weight.
    overlaps = mf.mo_coeff[:, :7].T @ mf.get_ovlp() @ mo_ordered[:, :10]
    from pyscf.fci import cistring
    determinants = cistring.make_strings(range(n), 3)
    slater_overlaps = np.array([
        np.linalg.det(overlaps[:, list(range(4)) + [4+i for i in range(n) if int(s) >> i & 1]])
        for s in determinants
    ])
    rhf_overlap = abs(slater_overlaps @ ci @ slater_overlaps)**2
    integral_path = args.folder / "integrals.npz"
    np.savez_compressed(integral_path, h1=h1, eri=eri, core=core,
                        fci_vector=vector, dm1=dm1, dm2=dm2,
                        mo_coeff=mo_ordered, ao_overlap=mf.get_ovlp())
    return {"distance_angstrom": args.distance, "molecule": "N2", "basis": "cc-pVDZ",
            "electrons": 14, "basis_orbitals": mol.nao, "active_electrons": 6,
            "active_orbitals": 6, "inactive_electrons": 8, "qubits": 12,
            "orbital_selection": "AVAS N 2p, ncore=4, threshold=0.2; singlet CASSCF(6,6)",
            "casscf_newton_fallback": newton_fallback,
            "active_orbital_permutation": order.tolist(),
            "rhf_energy": mf.e_tot, "casscf_energy": mc.e_tot, "fci_energy": efci,
            "core_energy": float(core), "rhf_ground_state_probability": rhf_overlap,
            "active_reference_determinant_probability": float(determinant_probability),
            "natural_occupations": np.linalg.eigvalsh(dm1)[::-1].tolist(),
            "spin_squared": float(solver.spin_square(ci, n, (3, 3))[0]),
            "integrals_sha256": hashlib.sha256(integral_path.read_bytes()).hexdigest()}


def hamiltonian(data):
    from qdk_chemistry.data import CanonicalFourCenterHamiltonianContainer, Hamiltonian, ModelOrbitals
    # The inactive orbitals are integrated into h1 and core. The quantum
    # register represents the effective active problem, not all 28 orbitals.
    return Hamiltonian(CanonicalFourCenterHamiltonianContainer(
        data["h1"], data["eri"].ravel(), ModelOrbitals(6), float(data["core"]), data["h1"]))


def contract(tensors, alpha_first):
    """Independent fermionic permutation, deliberately not adapter helpers."""
    n = len(tensors)
    result = np.zeros(4**n, dtype=complex)
    for local_states in itertools.product(range(4), repeat=n):
        value, modes = np.ones((1, 1)), []
        for site, (tensor, local) in enumerate(zip(tensors, local_states, strict=True)):
            value = value @ tensor[:, local, :]
            for spin in ((0, 1) if alpha_first else (1, 0)):
                if local & (1 << spin):
                    modes.append(spin*n + site)
        swaps = sum(x > y for i, x in enumerate(modes) for y in modes[i+1:])
        result[sum(1 << m for m in modes)] = (-1)**swaps * value[0, 0]
    return result


def solve(args):
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import Wavefunction
    import qdk_chemistry_alps  # noqa: F401

    data = np.load(args.folder / "integrals.npz")
    ham = hamiltonian(data)
    cas = create("multi_configuration_calculator", "macis_cas")
    cas.settings().set("ci_residual_tolerance", 1e-11)
    start = time.perf_counter()
    ecas, _ = cas.run(ham, 3, 3)
    cas_seconds = time.perf_counter() - start
    reference = json.loads((args.folder / "prepare.json").read_text())
    np.testing.assert_allclose(ecas, reference["fci_energy"], atol=1e-8, rtol=0)
    records = []
    for bond in args.bonds:
        calc = create("multi_configuration_calculator", "alps_classic_dmrg")
        for key, value in {"wavefunction_container": "mps", "max_bond_dim": bond,
                           "sweeps": args.sweeps,  "integral_cutoff": 1e-12,
                           "truncation_tol": 0.0, "convergence_tol": 1e-11}.items():
            calc.settings().set(key, value)
        start = time.perf_counter()
        energy, wfn = calc.run(ham, 3, 3)
        seconds = time.perf_counter() - start
        container = wfn.get_container()
        dm1 = np.asarray(wfn.get_active_one_rdm_spin_traced())
        dm2 = np.asarray(wfn.get_active_two_rdm_spin_traced()).reshape((6,)*4)
        erdm = float(np.einsum("pq,pq", data["h1"], dm1)
                     + .5*np.einsum("pqrs,pqrs", data["eri"], dm2) + data["core"])
        native = contract(calc.last_state.tensors, False)
        exported = contract([s.to_dense().reshape(s.shape) for s in container.sites], True)
        np.testing.assert_allclose(native, exported, atol=1e-11, rtol=0)
        np.testing.assert_allclose(erdm, energy, atol=1e-8, rtol=0)
        np.testing.assert_allclose(np.vdot(native, native), 1, atol=1e-9)
        from pyscf.fci.spin_op import spin_square0
        np.testing.assert_allclose(native.imag, 0, atol=1e-12)
        spin_squared = float(spin_square0(native[sector_indices()].real.reshape(20, 20), 6, (3, 3))[0])
        # PR #540 requires declared center zero. Verify the actual gauge
        # provided by the bridge; no state modification or relabeling
        # of noncanonical tensors is allowed.
        iso_error = max(np.linalg.norm((v := s.to_dense().reshape(s.shape)).reshape(s.shape[0], -1)
                                      @ v.reshape(s.shape[0], -1).T - np.eye(s.shape[0]))
                        for s in container.sites[1:])
        assert iso_error < 1e-9, iso_error
        canonical = wfn  # The bridge already supplies center zero and all spin RDMs.
        path = args.folder / f"d{bond}.wavefunction.h5"
        canonical.to_file(str(path), "hdf5")
        reloaded = Wavefunction.from_file(str(path), "hdf5")
        assert reloaded.content_hash() == canonical.content_hash()
        np.testing.assert_allclose(contract([s.to_dense().reshape(s.shape)
                                            for s in reloaded.get_container().sites], True), native, atol=1e-11)
        history = np.asarray([it.energy for it in calc.last_result.iterations])
        record = {"bond_requested": bond, "bond_actual": container.max_bond_dimension,
                  "energy": energy, "error_vs_exact_mEh": 1000*(energy-ecas), "seconds": seconds,
                  "rdm_energy_error": abs(erdm-energy), "one_rdm_max_error": float(abs(dm1-data["dm1"]).max()),
                  "fidelity_vs_exact": float(abs(np.vdot(data["fci_vector"], native))**2),
                  "spin_squared": spin_squared,
                  "natural_occupations": np.linalg.eigvalsh(dm1)[::-1].tolist(),
                  "last_iteration_delta_Eh": float(abs(history[-1]-history[-2])),
                  "iteration_energies": history.tolist(),
                  "convergence_note": "adjacent iterations can have different center positions; not a sweep convergence test",
                  "right_isometry_error": float(iso_error), "hdf5_bytes": path.stat().st_size,
                  "mps_tensor_bytes": sum(s.to_dense().nbytes for s in container.sites),
                  "serialization_hash": canonical.content_hash(), "process_peak_rss_bytes": rss()}
        records.append(record)
        print(json.dumps(record), flush=True)
    return {"casci_energy": ecas, "casci_seconds": cas_seconds, "dmrg": records}


def block2(args):
    from pyblock2.driver.core import DMRGDriver, SymmetryTypes
    data = np.load(args.folder / "integrals.npz")
    with tempfile.TemporaryDirectory(prefix="n2-block2-") as scratch:
        driver = DMRGDriver(scratch=scratch, symm_type=SymmetryTypes.SZ,
                            n_threads=1, stack_mem=500_000_000)
        driver.initialize_system(n_sites=6, n_elec=6, spin=0)
        start = time.perf_counter()
        mpo = driver.get_qc_mpo(h1e=data["h1"], g2e=data["eri"], ecore=float(data["core"]), iprint=0)
        ket = driver.get_random_mps(tag="N2", bond_dim=64, nroots=1)
        energy = driver.dmrg(mpo, ket, n_sweeps=24, bond_dims=[64],
                             noises=[1e-5]*4 + [1e-7]*4 + [0], thrds=[1e-12], tol=1e-11, iprint=0)
        dm1 = np.asarray(driver.get_1pdm(ket, iprint=0)).sum(axis=0)
        seconds = time.perf_counter()-start
        reference = json.loads((args.folder / "prepare.json").read_text())
        assert abs(energy-reference["fci_energy"]) < 1e-8
        assert abs(dm1-data["dm1"]).max() < 1e-6
        return {"energy": energy, "error_vs_exact_mEh": 1000*(energy-reference["fci_energy"]),
                "seconds_including_one_rdm": seconds, "symmetry": "Sz", "bond": 64,
                "one_rdm_max_error": float(abs(dm1-data["dm1"]).max()), "stack_bytes": 500_000_000}


def circuit(args):
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import MajoranaMapping, Wavefunction
    from qdk_chemistry.utils.qsharp import QSHARP_UTILS, create_qsharp_context, use_qsharp_context

    data = np.load(args.folder / "integrals.npz")
    wfn = Wavefunction.from_file(str(args.folder / f"d{args.bond}.wavefunction.h5"), "hdf5")
    target = contract([s.to_dense().reshape(s.shape) for s in wfn.get_container().sites], True)
    context = create_qsharp_context()
    with use_qsharp_context(context):
        context.set_quantum_seed(7)
        prep = create("state_prep", "matrix_product_state", unitary_synthesis=args.method,
                      rotation_bit_precision=args.bits)
        start = time.perf_counter()
        result = prep.run(wfn)
        build_seconds = time.perf_counter()-start
        # Persist the public synthesis parameters, sufficient to recreate the
        # operation under the pinned PR; avoid the PR's resetting QIR wrapper.
        parameters = prep.generate_matrix_product_state_preparation_data(wfn.get_container()).to_qsharp_params(args.bits)
        save(args.folder / f"circuit-d{args.bond}-{args.method}-{args.bits}bit-params.json", parameters)
        # Keep the virtual bond register until the dump. Tracing it out gives
        # the exact ensemble produced by the consumer's ResetAll(ancilla),
        # instead of silently selecting its most likely measurement outcome.
        namespace, opname, paramtype = (
            ("MPSSequential", "MPSSequential", "MPSSequentialParams") if args.method == "dense"
            else ("MPSSparse", "MPSSparse", "MPSSparseParams"))
        context.eval(f"""
        operation InspectNitrogenMPS(p : QDKChemistry.Utils.{namespace}.{paramtype}) : Unit {{
            use state = Qubit[p.numSites * p.numQubitsPerSite];
            use bond = Qubit[p.numAncillaQubits];
            {{
                use gradient = Qubit[p.rotationBitPrecision];
                within {{ QDKChemistry.Utils.PhaseGradient.PreparePhaseGradientState(gradient); }}
                apply {{ QDKChemistry.Utils.{namespace}.{opname}(
                    p.initialStateVec, p.numSites, p.siteToOrbitalOrder, p.siteDecompositions,
                    state, bond, gradient); }}
            }}
            Std.Diagnostics.DumpMachine();
            ResetAll(state + bond);
        }}
        """)
        typed_params = getattr(getattr(QSHARP_UTILS, namespace), paramtype)(**parameters)
        print(f"Circuit constructed: D={args.bond}, {args.bits} bits, {args.method}; simulating all bond outcomes", flush=True)
        start = time.perf_counter()
        simulation = context.run(context.code.InspectNitrogenMPS, 1, typed_params, save_events=True)
        dump = simulation[0]["events"][-1].state_dump().get_dict()
        bond_dim = 1 << parameters["numAncillaQubits"]
        prepared = np.zeros((1 << result.num_qubits, bond_dim), dtype=complex)
        for index, amplitude in dump.items():
            assert index < prepared.size, "Unexpected unreleased simulator ancilla"
            prepared[index // bond_dim, index % bond_dim] = amplitude
        prepared = prepared.reshape((2,)*result.num_qubits + (bond_dim,))
        prepared = prepared.transpose(*range(result.num_qubits-1, -1, -1), result.num_qubits).reshape(-1, bond_dim)
        simulate_seconds = time.perf_counter()-start
        np.savez_compressed(args.folder / f"circuit-d{args.bond}-{args.method}-{args.bits}bit-state.npz",
                            state_and_bond=prepared, simulate_seconds=simulate_seconds)
        print(f"Simulation saved after {simulate_seconds:.2f}s; mapping Hamiltonian and counting resources", flush=True)
        factory = result._qsharp_factory
        logical = context.logical_counts(factory.program, *factory.parameter.values())
        counts = dict(logical)
    np.testing.assert_allclose(np.sum(abs(prepared)**2), 1, atol=1e-9)
    # The QDK mapper returns the active electronic operator. Restore the
    # frozen/inactive+nuclear constant explicitly for total molecular energy.
    qham = create("qubit_mapper").run(hamiltonian(data), MajoranaMapping.jordan_wigner(12))
    from scipy.sparse import eye
    matrix = qham.to_matrix(sparse=True) + float(data["core"])*eye(4096, format="csr")
    exact = data["fci_vector"]
    reference_energy = json.loads((args.folder / "prepare.json").read_text())["fci_energy"]
    # QDK matrices use little-endian blocked modes; assert using an external
    # eigenvector before trusting the prepared circuit's energy.
    np.testing.assert_allclose(matrix @ exact, reference_energy*exact, atol=1e-8, rtol=0)
    energy = float(np.sum(prepared.conj() * (matrix @ prepared)).real)
    sector_probability = float(np.sum(abs(prepared[sector_indices()])**2))
    clean_bond_probability = float(np.sum(abs(prepared[:, 0])**2))
    return {"bond": args.bond, "method": args.method, "rotation_bits": args.bits,
            "state_qubits": result.num_qubits, "build_seconds": build_seconds,
            "simulate_seconds": simulate_seconds, "logical_counts": counts,
            "fidelity_vs_mps": float(np.sum(abs(target.conj() @ prepared)**2)),
            "fidelity_vs_exact": float(np.sum(abs(exact.conj() @ prepared)**2)),
            "zero_bond_probability": clean_bond_probability,
            "postselected_fidelity_vs_exact": float(abs(np.vdot(exact, prepared[:, 0]))**2 / clean_bond_probability),
            "energy": energy, "energy_error_mEh": 1000*(energy-reference_energy),
            "correct_particle_sector_probability": sector_probability,
            "pauli_terms": len(qham.pauli_strings), "wavefunction_hash": wfn.content_hash(),
            "execution": "Q# sparse simulator, MPSSequential/MPSSparse; all bond outcomes included, no postselection or hardware",
            "export_limit": "PR #540 circuit/QIR factory resets output; saved callable parameters used instead"}


def qpe_resources(args):
    """Compose the saved MPS with QDK QPE and count the resulting circuit.

    This is an explicit resource scenario, not a phase-estimation measurement
    or a certificate that the chosen Trotter error meets chemical accuracy.
    """
    from qdk_chemistry.algorithms import create
    from qdk_chemistry.data import AlgorithmRef, MajoranaMapping, QubitOperator, Wavefunction
    from qdk_chemistry.utils.qsharp import create_qsharp_context, use_qsharp_context
    data = np.load(args.folder / "integrals.npz")
    wfn = Wavefunction.from_file(str(args.folder / f"d{args.bond}.wavefunction.h5"), "hdf5")
    qham = create("qubit_mapper").run(hamiltonian(data), MajoranaMapping.jordan_wigner(12))
    # Shift using the solver's variational energy, not an exact CI answer.
    solved = json.loads((args.folder / "solve.json").read_text())
    total_guess = next(r["energy"] for r in solved["dmrg"] if r["bond_requested"] == args.bond)
    shift = total_guess - float(data["core"])
    coefficients = qham.coefficients.copy()
    coefficients[list(qham.pauli_strings).index("I"*12)] -= shift
    qham = QubitOperator(list(qham.pauli_strings), coefficients, encoding="jordan-wigner",
                        fermion_mode_order="blocked")
    context = create_qsharp_context()
    with use_qsharp_context(context):
        prep = create("state_prep", "matrix_product_state", unitary_synthesis=args.method,
                      rotation_bit_precision=args.bits).run(wfn)
        time_scale, phase_bits, divisions = 0.5, args.phase_bits, args.trotter_divisions
        builder = create("qpe_circuit_builder", "qdk_standard", num_bits=phase_bits,
                         controlled_circuit_mapper=AlgorithmRef("controlled_circuit_mapper", "pauli_sequence"),
                         unitary_builder=AlgorithmRef("hamiltonian_unitary_builder", "trotter",
                                                      time=time_scale, order=2, num_divisions=divisions))
        full = builder.run(state_preparation=prep, qubit_hamiltonian=qham)[0]
        factory = full._qsharp_factory
        counts = dict(context.logical_counts(factory.program, *factory.parameter.values()))
    return {"bond": args.bond, "rotation_bits": args.bits, "phase_bits": phase_bits,
            "evolution_time": time_scale, "trotter_order": 2, "trotter_divisions": divisions,
            "energy_shift_electronic_Eh": shift, "core_energy_Eh": float(data["core"]),
            "energy_bin_width_Eh": 2*np.pi/(time_scale*2**phase_bits), "logical_counts": counts,
            "wavefunction_hash": wfn.content_hash(),
            "status": "Composed and resource-counted in QDK; QPE not executed",
            "limitation": "Phase bin width is not total error. Trotter and gate-synthesis errors are not certified; no physical-qubit or runtime claim."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "solve", "block2", "circuit", "qpe_resources"])
    parser.add_argument("--folder", type=Path, required=True)
    parser.add_argument("--distance", type=float, default=1.1)
    parser.add_argument("--bonds", type=int, nargs="+", default=[8, 16, 32, 64])
    parser.add_argument("--sweeps", type=int, default=12)
    parser.add_argument("--bond", type=int, default=16)
    parser.add_argument("--bits", type=int, default=8)
    parser.add_argument("--method", choices=["dense", "block_sparse"], default="dense")
    parser.add_argument("--phase-bits", type=int, default=8, help="QPE resolution bits; default is a small resource-counting demo")
    parser.add_argument("--trotter-divisions", type=int, default=4, help="QPE second-order Trotter subdivisions")
    args = parser.parse_args()
    if args.phase_bits < 1 or args.trotter_divisions < 1:
        parser.error("phase-bits and trotter-divisions must be positive")
    args.folder.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    record = globals()[args.stage](args)
    record.update(total_seconds=time.perf_counter()-start, process_peak_rss_bytes=rss(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  python=sys.version, versions={p: importlib.metadata.version(p)
                                               for p in ("numpy", "scipy", "pyscf", "qdk")})
    suffix = f"-d{args.bond}-{args.method}-{args.bits}bit" if args.stage in ("circuit", "qpe_resources") else ""
    save(args.folder / f"{args.stage}{suffix}.json", record)
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
