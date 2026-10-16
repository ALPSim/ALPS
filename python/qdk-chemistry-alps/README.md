# qdk-chemistry-alps

ALPS solver backends for Microsoft's [qdk-chemistry](https://github.com/microsoft/qdk-chemistry), built on [pyalps](../pyalps/README.md): exact diagonalization and classic block DMRG registered as `MultiConfigurationCalculator`s, plus adapters between the two packages' lattice, Hamiltonian and wavefunction representations.

```python
import qdk_chemistry_alps                    # registers the backends
from qdk_chemistry.algorithms import create

energy, wavefunction = create("multi_configuration_calculator", "alps_ed").run(hamiltonian, n_alpha, n_beta)

dmrg = create("multi_configuration_calculator", "alps_classic_dmrg", max_bond_dim=200, sweeps=12)
energy, wavefunction = dmrg.run(hamiltonian, n_alpha, n_beta)
dm1, dm2 = dmrg.make_rdms()                  # spin-summed, PySCF convention
```

## What it provides

| | |
|---|---|
| **`alps_ed`** | Exact diagonalization of the integrals as a `pyalps.models` Hamiltonian, one spinful-fermion site per orbital, in the (N, S<sub>z</sub>) sector with `pyalps.ed`. Returns the CI vector as a `StateVectorContainer`. Settings: `integral_cutoff` (1e-12), `ci_coefficient_cutoff` (0: keep every amplitude; a positive value drops smaller ones and renormalizes). |
| **`alps_classic_dmrg`** | Classic block DMRG of restricted real integrals through `pyalps.classic_dmrg.solve_molecular`, which builds the Hamiltonian, sweeps and measures every 1-/2-RDM element in C++. Returns an `RdmContainer` wavefunction (spin-resolved 1-RDMs, spin-summed 2-RDM); `make_rdms()` and `last_result` (iterations, truncation error) stay on the calculator. Registered when pyalps has its `classic_dmrg` bindings. |
| **`interop`** | `lattice_graph_to_alps` / `alps_lattice_to_lattice_graph`, `qubit_operator_to_hamiltonian` / `hamiltonian_to_qubit_operator` (spin-1/2), `pauli_operator_to_qubit_operator` (the native `pyalps.models.export_pauli` route: finite spins and bosons with binary or Gray encoding, strict fermions with Jordan–Wigner strings, pairing included), `hamiltonian_to_integrals` and `fermionic_hamiltonian_to_qubit_operator` (restricted fermionic models through QDK's Jordan–Wigner, Bravyi–Kitaev, BK-tree and parity mappings), `ci_vector_to_qubit_statevector`, `rcm_ordering`. |

### `alps_classic_dmrg` settings

| Setting | Default | Meaning |
| --- | ---: | --- |
| `max_bond_dim` | 200 | Kept states in every half sweep, and the cap under error-based truncation. |
| `warmup_states` | 20 | Kept states while the blocks grow. |
| `sweeps` | 20 | Sweeps before the final measurement pass. |
| `convergence_tol` | 1e-10 | Lanczos tolerance of the engine; inspect `last_result.iterations` for sweep convergence. |
| `truncation_tol` | 0 | Zero keeps the fixed-state schedule; a positive value targets that discarded weight, capped at `max_bond_dim`. |
| `integral_cutoff` | 1e-12 | Integrals of smaller magnitude are dropped before the symmetry check. |

Fixed-state sweeps are the default: alpsx reports that in its LiH validation error-based truncation stopped at a higher energy despite a tiny reported discarded weight, so check energy and RDM convergence when opting in. The solver suits small active spaces. Every integral is a separate engine term, regrouped at every step rather than kept in complementary operators, and the RDM measurement adds n⁴ more terms. With 32 kept states and two sweeps on dense integrals, an Apple M-series laptop takes about 1.5 s for 6 orbitals, 9 s for 8, 30 s for 10 and 100 s for 12 (see the [component README](../../src/alps/classic_dmrg/README.md#molecular-hamiltonians-alpsdmrgmolecularhpp)). The result holds no MPS and no determinant coefficients, so there are no later measurements, no CI reconstruction and no continuation.

## Conventions

- **Mode order.** ALPS orders a site's modes down before up (`|up,dn> = c+_dn c+_up |0>`) and the sites ascending; qdk-chemistry, MACIS and PySCF order all alpha modes before all beta modes. `abinitio.interleaved_to_blocked_parity` is the sign between the two, applied to every CI coefficient `alps_ed` returns; the tests pin it against PySCF's FCI vector.
- **Indices.** QDK Pauli strings put qubit 0 at the right end and its dense matrices put qubit 0 in the least significant bit; pyalps product bases put site 0 in the most significant position, and its Pauli export puts site 0 in the least significant register bits. `ci_vector_to_qubit_statevector` places a CI vector in the blocked Jordan–Wigner register of the integral route.
- **RDMs.** `dm1[p, q] = <c+_q c_p>` and `dm2[p, q, r, s] = sum_ab <c+_pa c+_rb c_sb c_qa>`, as in PySCF and QDK's `RdmContainer`; `E = sum h dm1 + 1/2 sum (pq|rs) dm2 + E_core`.

## Installing

pyalps is not a declared dependency: the backends need its `models`, `ed` and `classic_dmrg` bindings, which a pyalps built from this repository against an installed ALPS SDK provides. Build the SDK (with `ALPS_ENABLE_CLASSIC_DMRG`, the default, for `alps_classic_dmrg`), install pyalps against it as its [README](../pyalps/README.md) describes, then

```bash
python -m pip install -e python/qdk-chemistry-alps
```

Without pyalps the import fails with this guidance. `alps_classic_dmrg` needs a qdk-chemistry build providing `RdmContainer`, which alpsx's compatibility notes list as missing from released versions up to 2.2.1, and raises `RuntimeError` without it; `alps_ed` and `interop` do not use it. This port was tested only with a local qdk-chemistry 2.1.0 build that has `RdmContainer`. Import registers the backends; qdk-chemistry builds with plugin discovery also load the `qdk_chemistry.plugins` entry point, and both paths are idempotent.

## Tests

```bash
python -m pytest python/qdk-chemistry-alps --require-solver-stack
```

The suite takes seconds. It compares `alps_ed` with MACIS (energies and signed CI vectors) and PySCF FCI (H₂, H₄, the sign-exact CI vector), `alps_classic_dmrg` with `alps_ed`, MACIS and PySCF FCI (energies and every 1-/2-RDM element, LiH in STO-3G with the default schedule), the integral extraction with the builder and with composite operators, the qubit export with the Fock-space spectrum for every encoding and with QDK's own mapper term by term, the state placement with the energy (and a bit-reversed register as a negative control), and the native Pauli route with independently built encodings. Markers gate `pyscf`, `local_qdk` (`RdmContainer`) and `classic_dmrg`; `--require-solver-stack` fails instead of skipping when pyalps, qdk-chemistry or the classic bindings are missing.

## Origin

Adapted from the `qdk-chemistry-alps` package of the private skilledwolf/alpsx fork (branch `agent/v3-execution-and-api-audit` at d37e51227, and the `alps_classic_dmrg` calculator of branch `WIP/classic-dmrg-chemistry` at d4384a7f8), which ran on that fork's unpublished `alps` Python package. Every file names its source and changes. Not carried over: the MPS DMRG calculator `alps_dmrg` and its measurements and MPS/state-vector containers (ALPS 3 has no MPS library), the documentation site, examples, figures and proof-of-concept scripts.

## License

MIT — see [LICENSE](LICENSE). To cite ALPS, see the repository's [CITATION.cff](../../CITATION.cff).
