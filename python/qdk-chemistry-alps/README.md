# qdk-chemistry-alps

ALPS solver backends for Microsoft's [qdk-chemistry](https://github.com/microsoft/qdk-chemistry), built on [pyalps](../pyalps/README.md): exact diagonalization and classic block DMRG registered as `MultiConfigurationCalculator`s, plus adapters between the two packages' lattice, Hamiltonian and wavefunction representations.

```python
import qdk_chemistry_alps                    # registers the backends
from qdk_chemistry.algorithms import create

energy, wavefunction = create("multi_configuration_calculator", "alps_ed").run(hamiltonian, n_alpha, n_beta)

dmrg = create("multi_configuration_calculator", "alps_classic_dmrg", max_bond_dim=200, sweeps=12)
energy, wavefunction = dmrg.run(hamiltonian, n_alpha, n_beta)   # CI vector and RDMs
dm1, dm2 = dmrg.make_rdms()                  # spin-summed, PySCF convention
spin = dmrg.spin_spin_correlations()         # <S_i . S_j>, from the RDMs
```

## What it provides

| | |
|---|---|
| **`alps_ed`** | Exact diagonalization of the integrals as a `pyalps.models` Hamiltonian, one spinful-fermion site per orbital, in the (N, S<sub>z</sub>) sector with `pyalps.ed`. Returns the CI vector as a `StateVectorContainer`. Settings: `integral_cutoff` (1e-12), `ci_coefficient_cutoff` (0: keep every amplitude; a positive value drops smaller ones and renormalizes). |
| **`alps_classic_dmrg`** | Classic block DMRG of restricted real integrals through `pyalps.classic_dmrg.solve_molecular`, which builds the Hamiltonian, sweeps and reads in C++, from the final state, every 1-/2-RDM element by spin block and, for small active spaces, the CI vector, contracted from the engine's stored rotations. Returns a `StateVectorContainer` with the CI vector and the run's RDMs, or an `RdmContainer` (see `wavefunction_container`). `make_rdms()`, `measure_local()`, `measure_correlations()`, `spin_spin_correlations()` and `last_result` (iterations, truncation error, observables) stay on the calculator. Registered when pyalps has its `classic_dmrg` bindings. |
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
| `wavefunction_container` | `auto` | `auto`: a `StateVectorContainer` while the CI vector fits in `ci_memory_limit_gb`, otherwise an `RdmContainer`. `statevector`: always the CI vector (`MemoryError` beyond the limit). `rdm`: an `RdmContainer` without a CI vector. |
| `ci_memory_limit_gb` | 4 | Budget for the CI vector, estimated at 1 kB per determinant of the sector (about 460 B measured for QDK's objects, plus the engine's arrays), so about four million determinants. The choice is made before the solve. |
| `ci_coefficient_cutoff` | 0 | Zero keeps every determinant; a positive value drops smaller amplitudes and renormalizes, and the container then carries no RDMs, as the trimmed vector's differ from the run's. |

Fixed-state sweeps are the default: alpsx reports that in its LiH validation error-based truncation stopped at a higher energy despite a tiny reported discarded weight, so check energy and RDM convergence when opting in. The engine keeps the normal and complementary block operators of quantum-chemistry DMRG and reads the RDMs by a walk over the final state (see the [component README](../../src/alps/classic_dmrg/README.md#molecular-hamiltonians-alpsdmrgmolecularhpp)). With 32 kept states and two sweeps on dense random integrals at half filling, an Apple M-series laptop takes 0.6 s for 8 orbitals, 1.4 s for 10 and 2.5 s for 12 in `pyalps.classic_dmrg.solve_molecular`. The CI vector is what limits `auto`: at 12 orbitals (853776 determinants) the engine reconstructs it in 0.35 s and the `StateVectorContainer` takes 2.4 s to build, 2.1 s of it in qdk-chemistry's own `Configuration` objects and container, so the calculator takes 5.8 s there against 2.9 s with `wavefunction_container="rdm"`.

The CI vector is the DMRG state, truncated or not: its energy and RDMs are the run's. The engine keeps no state after the run, so everything is measured in the final sweep. `measure_local` (`n`, `n_up`, `n_dn`, `Sz`, `n_updn`), `measure_correlations` (`<O_i+ O_j>` for the number-type operators `n`, `n_up`, `n_dn` and `Sz` in any combination, or one of `c_up`, `c_dn`, `Sp` and `Sm` on both sides) and `spin_spin_correlations` derive their values from the RDMs, as alpsx's `alps_dmrg` contracted them from its MPS. For other operators, set `calculator.observables` to a dict of `pyalps.models` term sums on `spinful_fermion()` sites (site `p` is orbital `p`, at most four sites per product) before `run()` and read `last_result.observables`. There is no continuation and no MPS export.

## Conventions

- **Mode order.** ALPS orders a site's modes down before up (`|up,dn> = c+_dn c+_up |0>`) and the sites ascending; qdk-chemistry, MACIS and PySCF order all alpha modes before all beta modes. `abinitio.interleaved_to_blocked_parity` is the sign between the two, applied to every CI coefficient `alps_ed` and `alps_classic_dmrg` return; the tests pin it against PySCF's FCI vector and MACIS.
- **Indices.** QDK Pauli strings put qubit 0 at the right end and its dense matrices put qubit 0 in the least significant bit; pyalps product bases put site 0 in the most significant position, and its Pauli export puts site 0 in the least significant register bits. `ci_vector_to_qubit_statevector` places a CI vector in the blocked Jordan–Wigner register of the integral route.
- **RDMs.** `dm1[p, q] = <c+_q c_p>` and `dm2[p, q, r, s] = sum_ab <c+_pa c+_rb c_sb c_qa>`, as in PySCF and QDK's `RdmContainer`; `E = sum h dm1 + 1/2 sum (pq|rs) dm2 + E_core`. QDK's `aaaa`, `aabb` and `bbbb` blocks are PySCF's `make_rdm12s` blocks, `aabb[p, q, r, s] = <c+_pa c+_rb c_sb c_qa>` with a = alpha, b = beta (checked against MACIS).

## Installing

pyalps is not a declared dependency: the backends need its `models`, `ed` and `classic_dmrg` bindings, which a pyalps built from this repository against an installed ALPS SDK provides. Build the SDK (with `ALPS_ENABLE_CLASSIC_DMRG`, the default, for `alps_classic_dmrg`), install pyalps against it as its [README](../pyalps/README.md) describes, then

```bash
python -m pip install -e python/qdk-chemistry-alps
```

Without pyalps the import fails with this guidance. `alps_classic_dmrg`'s `rdm` container needs a qdk-chemistry build providing `RdmContainer`, which the released 2.2.1 lacks (alpsx's compatibility notes report the same for earlier releases); there it returns state vectors, and raises `RuntimeError` naming the class before solving when the settings call for an `RdmContainer`. `alps_ed` and `interop` do not use it. The suite passes with qdk-chemistry 2.2.1 from PyPI on Python 3.13, where the tests of the `RdmContainer` path skip, and in full with a local 2.1.0 build that has `RdmContainer`. Import registers the backends; qdk-chemistry builds with plugin discovery, 2.2.1 among them, also load the `qdk_chemistry.plugins` entry point, and both paths are idempotent.

## Tests

```bash
python -m pytest python/qdk-chemistry-alps --require-solver-stack
```

The suite takes seconds. It compares `alps_ed` with MACIS (energies and signed CI vectors) and PySCF FCI (H₂, H₄, the sign-exact CI vector), `alps_classic_dmrg` with `alps_ed`, MACIS and PySCF FCI (energies, signed CI vectors, every 1-/2-RDM element and spin block, LiH in STO-3G with the default schedule; a truncated run's CI vector against its own energy; the measurements and an observable against `pyalps.ed`, and ⟨S²⟩ against QDK's), the integral extraction with the builder and with composite operators, the qubit export with the Fock-space spectrum for every encoding and with QDK's own mapper term by term, the state placement with the energy (and a bit-reversed register as a negative control), and the native Pauli route with independently built encodings. Markers gate `pyscf`, `local_qdk` (`RdmContainer`) and `classic_dmrg`; `--require-solver-stack` fails instead of skipping when pyalps, qdk-chemistry or the classic bindings are missing.

## Origin

Adapted from the `qdk-chemistry-alps` package of the private skilledwolf/alpsx fork (branch `agent/v3-execution-and-api-audit` at d37e51227, and the `alps_classic_dmrg` calculator of branch `WIP/classic-dmrg-chemistry` at d4384a7f8), which ran on that fork's unpublished `alps` Python package. Every file names its source and changes. `alps_classic_dmrg`'s container modes, CI settings and measurement methods follow `alps_dmrg` (`dmrg_calculator.py` at d37e51227), with the CI vector and the measurements taken from the classic engine. Not carried over: the MPS DMRG calculator `alps_dmrg` itself and its MPS container (ALPS 3 has no MPS library), the documentation site, examples, figures and proof-of-concept scripts.

## License

MIT — see [LICENSE](LICENSE). To cite ALPS, see the repository's [CITATION.cff](../../CITATION.cff).
