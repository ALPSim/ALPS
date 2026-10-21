# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/spin_one_quench.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 4d00f034c76644883a29a0dac559b00fc8df0f7f60f48be9f864788fed01c0dc.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Spin-1 quench: physical transport, executed Q# circuits and encoding leakage."""

import time
from collections import Counter

import numpy as np
from demo_support import export_operator
from pyalps import models
from demo_support import model_matrix
from pyalps import classic_dmrg as dmrg
from pyalps.models import make_qubit_encoding
from demo_support import bit_reversal, finish, parser
from demo_support.spin_one import (
    COUPLINGS,
    FIELDS,
    model,
    qubit_matrix,
    reference_indices,
)
from qdk import qsharp

from qdk_chemistry_alps.interop import pauli_operator_to_qubit_operator

TIMES = np.arange(5)*.2
INITIAL = [0, 1, 2, 1]  # m = +1, 0, -1, 0
STEPS = [.1, .05, .025]
TARGETS = {"maximum_observable_error": 5e-4, "maximum_leakage": 1e-7, "maximum_infidelity": 1e-6}


def source_for(encoded, encoding):
    gates = []
    weights = []
    for word, coefficient in encoded.terms.items():
        assert abs(coefficient.imag) < 1e-12
        if set(word) == {"I"}:
            continue  # irrelevant global phase for these observables/fidelities
        paulis = ", ".join("Pauli"+p for p in word[::-1])
        gates.append(f"Exp([{paulis}], ({-coefficient.real:.17f}) * dt / 2.0, q);")
        weights.append(len(word)-word.count("I"))
    prepare = " ".join(f"X(q[{q}]);" for q, b in enumerate(encoding.encode_basis_bits(INITIAL)) if b)
    source = """operation SpinOnePrepare(q : Qubit[]) : Unit is Adj + Ctl {
    """ + prepare + """
}
operation SpinOneEvolve(q : Qubit[], steps : Int, dt : Double) : Unit is Adj + Ctl {
    for _ in 1..steps {
        """ + "\n        ".join(gates + list(reversed(gates))) + """
    }
}
operation SpinOneShot(steps : Int, dt : Double) : Result[] {
    use q = Qubit[8];
    SpinOnePrepare(q);
    SpinOneEvolve(q, steps, dt);
    let bits = [M(q[0]), M(q[1]), M(q[2]), M(q[3]), M(q[4]), M(q[5]), M(q[6]), M(q[7])];
    ResetAll(q);
    return bits;
}
"""
    resources = {"pauli_rotations_per_step": 2*len(gates),
                 "cnot_ladder_estimate_per_step": 4*sum(w-1 for w in weights),
                 "initial_x_gates": sum(encoding.encode_basis_bits(INITIAL)),
                 "pauli_weight_histogram_without_identity": dict(sorted(Counter(weights).items())),
                 "resource_assumption": "Unoptimized all-to-all CNOT ladder for each Exp; not compiled hardware counts. Global phase omitted."}
    return source, resources


def trajectory(source, dt):
    qsharp.init()
    qsharp.eval(source)
    # StateDump order is independently probed; chemistry Pauli matrices use
    # q0 least significant, while this interpreter dumps q0 most significant.
    qsharp.eval("use q = Qubit[8]; X(q[0]);")
    probe = np.asarray(qsharp.dump_machine().as_dense_state())
    assert abs(probe[128]-1.) < 1e-13
    qsharp.eval("X(q[0]); SpinOnePrepare(q);")
    reverse = bit_reversal(8)
    states = []
    for index, _ in enumerate(TIMES):
        if index:
            qsharp.eval(f"SpinOneEvolve(q, {round(.2/dt)}, {dt:.17f});")
        states.append(np.asarray(qsharp.dump_machine().as_dense_state())[reverse])
    qsharp.eval("ResetAll(q);")
    return np.asarray(states)


def mps_trajectory(h, dt, bond):
    spin = h.site_for(0)
    def z(i):
        return models.located(spin, "Sz", i)
    result = dmrg.evolve(INITIAL, h, time_step=dt, steps=round(TIMES[-1]/dt),
                         states=bond, method="trotter",
                         measure_steps=list(range(0, round(TIMES[-1]/dt)+1, round(.2/dt))),
                         measurements={"Sz": sum((z(i) for i in range(4)), start=models.TermSum()),
                                       "Sz2": sum((z(i)*z(i) for i in range(4)), start=models.TermSum()), "Czz": z(0)*z(1)})
    np.testing.assert_allclose(result.times, TIMES, atol=1e-12)
    return np.column_stack([result.measurements[k].values.real for k in ("Sz", "Sz2", "Czz")])


def plot(plt, data):
    figure, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    exact = np.asarray(data["exact_observables"])
    for site, color in [(0, "#087e8b"), (2, "#7956a8")]:
        axes[0].plot(TIMES, exact[:, site], color=color, label=f"Site {site}: exact")
        axes[0].plot(TIMES, np.asarray(data["mps_fine"])[:, site], "+", color=color)
        for scheme, marker in [("binary", "o"), ("gray", "x")]:
            values = np.asarray(data["encodings"][scheme]["runs"][-1]["observables"])
            axes[0].plot(TIMES, values[:, site], marker, fillstyle="none", color=color, markersize=4)
    axes[0].plot([], [], "+", color="0.3", label="MPS")
    axes[0].plot([], [], "o", fillstyle="none", color="0.3", label="Q# binary")
    axes[0].plot([], [], "x", color="0.3", label="Q# Gray")
    axes[0].set(xlabel="Time (ℏ/J₀)", ylabel="⟨Sᶻ⟩", title="Spin transport: exact, MPS and Q#")
    axes[0].legend(fontsize=8)
    axes[1].plot(TIMES, exact[:, 5], color="#087e8b", label="Exact, site 1")
    for scheme, marker, color in [("binary", "o", "#d68c22"), ("gray", "x", "#7956a8")]:
        entry = data["encodings"][scheme]
        axes[1].plot(TIMES, np.asarray(entry["runs"][-1]["observables"])[:, 5], marker,
                     color=color, fillstyle="none", label=f"Q# {scheme}")
        axes[2].loglog(STEPS, [max(r["leakage"])+1e-16 for r in entry["runs"]], marker+"-", color=color, label=scheme)
    axes[1].set(xlabel="Time (ℏ/J₀)", ylabel="⟨(Sᶻ₁)²⟩", title="Exchange populates m = ±1")
    axes[1].legend(fontsize=8)
    axes[2].set(xlabel="Q# product-formula step Δt", ylabel="Maximum probability outside physical space",
                title="Smaller steps reduce leakage")
    axes[2].legend(fontsize=8)
    return figure


def main():
    arguments = parser(__doc__)
    arguments.add_argument("--shots", type=int, default=64, help="Q# shots at the final time per encoding")
    args = arguments.parse_args()
    if args.shots < 1:
        arguments.error("--shots must be positive")
    started = time.perf_counter()
    h, reference, axes = model(.65, biquadratic=0., transverse=.17, y_field=.09)
    np.testing.assert_allclose(model_matrix(h), reference, atol=1e-13)
    initial = np.zeros(81, dtype=complex)
    initial[np.ravel_multi_index(tuple(INITIAL), (3, 3, 3, 3))] = 1.
    energies, vectors = np.linalg.eigh(reference)
    exact_states = (np.exp(-1j*TIMES[:, None]*energies)*(vectors.conj().T @ initial)) @ vectors.T
    physical_observables = axes[2] + [z@z for z in axes[2]] + [axes[2][0]@axes[2][1]]
    exact = np.column_stack([np.einsum("ti,ij,tj->t", exact_states.conj(), o, exact_states).real for o in physical_observables])
    classical_coarse = mps_trajectory(h, .025, 9)
    classical_fine = mps_trajectory(h, .0125, 9)
    mps_errors = {"dt_0.025_bond_9": float(np.max(abs(classical_coarse-exact))),
                  "dt_0.0125_bond_9": float(np.max(abs(classical_fine-exact)))}
    assert mps_errors["dt_0.0125_bond_9"] < mps_errors["dt_0.025_bond_9"]*.4
    np.testing.assert_allclose(classical_fine, exact, atol=.001, rtol=0)
    spin = h.site_for(0)
    expressions = ([models.located(spin, "Sz", i) for i in range(4)] + [models.located(spin, "Sz", i)*models.located(spin, "Sz", i) for i in range(4)]
                   + [models.located(spin, "Sz", 0)*models.located(spin, "Sz", 1)])
    data = {"times": TIMES, "exact_observables": exact, "mps_fine": classical_fine,
            "mps_coarse": classical_coarse,
            "mps_errors": mps_errors, "encodings": {},
            "observable_order": [f"Sz[{i}]" for i in range(4)] + [f"Sz2[{i}]" for i in range(4)] + ["Sz[0] Sz[1]"]}
    print("Four-site spin-1 quench: actual Q# execution, both compact encodings, no postselection.")
    for number, scheme in enumerate(("binary", "gray")):
        encoding = make_qubit_encoding(h, scheme=scheme)
        encoded = export_operator(h, encoding=encoding, tolerance=0.)
        matrix = qubit_matrix(pauli_operator_to_qubit_operator(encoded))
        indices = reference_indices(scheme)
        V = np.eye(256)[:, indices]
        np.testing.assert_allclose(matrix @ V, V @ reference, atol=1e-12)
        observables = [qubit_matrix(pauli_operator_to_qubit_operator(
            export_operator(h, observable=o, encoding=encoding, tolerance=0.))) for o in expressions]
        for actual, expected in zip(observables, physical_observables, strict=True):
            np.testing.assert_allclose(actual @ V, V @ expected, atol=1e-12)
        # Sz and Sz² are diagonal even on unused states; use their actual
        # encoded values for unconditional shot estimates.
        diagonals = np.column_stack([o.diagonal().real for o in observables])
        valid = np.zeros(256, dtype=bool)
        valid[indices] = True
        mapped_exact = np.asarray([encoding.encode_statevector(state) for state in exact_states])
        source, resources = source_for(encoded, encoding)
        runs = []
        for dt in STEPS:
            states = trajectory(source, dt)
            np.testing.assert_allclose(np.linalg.norm(states, axis=1), 1., atol=1e-11)
            values = abs(states)**2 @ diagonals
            # Measure a fraction of the actual norm, without projecting or
            # renormalizing states. Avoid subtractive cancellation near zero.
            norm_squared = np.sum(abs(states)**2, axis=1)
            leakage = np.sum(abs(states[:, ~valid])**2, axis=1) / norm_squared
            np.testing.assert_allclose(leakage, [1-encoding.physical_weight(s) for s in states], atol=1e-12)
            fidelity = abs(np.sum(mapped_exact.conj()*states, axis=1))**2
            runs.append({"dt": dt, "observables": values, "leakage": leakage, "fidelity": fidelity,
                         "norm_squared": norm_squared,
                         "maximum_observable_error": float(np.max(abs(values-exact))),
                         "steps_to_final_time": round(TIMES[-1]/dt),
                         "cnot_ladder_estimate": resources["cnot_ladder_estimate_per_step"]*round(TIMES[-1]/dt)})
        assert runs[-1]["maximum_observable_error"] < .4*runs[-2]["maximum_observable_error"]
        assert max(runs[-1]["leakage"]) < max(runs[-2]["leakage"])*.4 + 1e-12
        eligible = [r for r in runs if r["maximum_observable_error"] < TARGETS["maximum_observable_error"]
                    and max(r["leakage"]) < TARGETS["maximum_leakage"]
                    and 1-min(r["fidelity"]) < TARGETS["maximum_infidelity"]]
        assert eligible, f"{scheme} did not reach the common accuracy targets"
        matched = min(eligible, key=lambda r: r["cnot_ladder_estimate"])
        seed = 971 + number
        rows = qsharp.run(f"SpinOneShot({round(TIMES[-1]/STEPS[-1])}, {STEPS[-1]})", shots=args.shots, seed=seed)
        if any(not isinstance(row, list) or len(row) != 8 for row in rows):
            raise RuntimeError(f"QDK shot failed: {rows!r}")
        measured = np.array([sum(int(bit == qsharp.Result.One) << q for q, bit in enumerate(row)) for row in rows])
        means = diagonals[measured].mean(axis=0)
        shot_bound = np.sqrt(2*np.log(2*len(means)/.001)/args.shots)
        assert np.max(abs(means-runs[-1]["observables"][-1])) <= shot_bound
        physical_fraction = valid[measured].mean()
        validity_bound = np.sqrt(np.log(2/.001)/(2*args.shots))
        assert abs(physical_fraction-(1-runs[-1]["leakage"][-1])) <= validity_bound
        data["encodings"][scheme] = {"layout": {"scheme": scheme, "site_dimensions": list(encoding.site_dimensions), "site_offsets": list(encoding.site_offsets)}, "pauli_terms": len(encoded.terms),
            "resources": resources, "runs": runs, "qsharp_source": source,
            "coefficient_l1": encoded.coefficient_l1(), "discarded_l1": encoded.discarded_l1,
            "matched_accuracy_dt": matched["dt"], "matched_accuracy_cnot_estimate": matched["cnot_ladder_estimate"],
            "matched_accuracy_pauli_rotations": resources["pauli_rotations_per_step"]*matched["steps_to_final_time"],
            "shot_count": args.shots, "seed": seed, "shot_observables": means,
            "shot_register_counts": {str(k): v for k, v in Counter(measured.tolist()).items()},
            "shot_physical_fraction": physical_fraction, "familywise_shot_bound": shot_bound,
            "shot_validity_bound": validity_bound, "pruning_tolerance": 0.}
        print(f"{scheme}: max observable error {runs[-1]['maximum_observable_error']:.3g}; "
              f"max leakage {max(runs[-1]['leakage']):.3g}; min fidelity {min(runs[-1]['fidelity']):.9f}")
        print(f"  Shared accuracy targets reached at dt={matched['dt']}: "
              f"{matched['cnot_ladder_estimate']} estimated CNOTs in the unoptimized ladder construction.")
    data["model"] = {"sites": 4, "spin": 1, "J": COUPLINGS, "K": 0., "D": .65, "h": FIELDS,
                     "x_field_site_1": .17, "y_field_site_2": .09, "constant": .37,
                     "initial_local_indices": INITIAL, "boundary": "open", "units": "hbar=1; time in hbar/J0"}
    data["controls"] = {"qdk_steps": STEPS, "qdk_order": 2, "mps_steps": [.025, .0125],
                        "mps_bonds": [9], "mps_order": 2, "truncation": None,
                        "bond_note": "Four sites remain explicit in the classic superblock; this example tests time-step error, not bond truncation",
                        "shots": args.shots, "postselection": False, "global_phase_omitted": True}
    data["accuracy_targets"] = TARGETS
    finish(args, data, started, plot)


if __name__ == "__main__":
    main()
