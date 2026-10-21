# Advanced QDK–ALPS examples

Start with [the reproducible installation guide](INSTALL.md). It selects the public
ALPS branch, builds the native bindings and pins the combined QDK MPS snapshot.
The adapter extra alone does not build pyalps or select that unreleased QDK API.

| Example | Calculation and independent check |
| --- | --- |
| `quickstart.py`, `molecule.py` | QDK calculator registration, Hubbard and LiH exact diagonalization |
| `interop.py`, `hardware_export.py` | Lattice round trips, Pauli operators and small spectra |
| `custom_spin_export.py` | Custom spin Hamiltonian and qubit conventions |
| `fermion_pairing.py` | Jordan–Wigner export with pairing; independent CAR reference |
| `hubbard_reference.py` | Molecular-integral and native fermionic routes against Fock space |
| `observables.py` | Local and two-site observables checked against stored RDM identities |
| `rdm_observables.py` | 14-orbital Hubbard RDMs and derived observables without CI reconstruction |
| `large_active_space.py` | 14-orbital QDK MPS export and energy from its stored RDMs |
| `converging_dmrg.py` | Kept-state scan against MACIS; final-state energies |
| `spin_one_anisotropy.py` | Four spin-1 sites, binary/Gray encodings and physical-subspace spectra |
| `spin_quench.py` | Six spin-1/2 sites: classic time evolution, exact dynamics and actual Q# execution |
| `spin_one_quench.py` | Spin-1 dynamics, encoding leakage and unconditional Q# shot estimates |
| [nitrogen](nitrogen/README.md) | Bond breaking, MPS serialization, Q# state preparation, QPE resource counts and CAS(10,26) capacity |
| [scaling](scaling/README.md) | Hubbard and full-Coulomb hydrogen chains with resource guards |

The detailed demonstrations (`custom_spin_export`, `fermion_pairing`,
`hubbard_reference`, `rdm_observables`, `spin_one_anisotropy`,
`spin_quench`, `spin_one_quench`) accept `--output results/demos` to save
JSON data, provenance and a figure. The other top-level examples print results. Run from
a normal terminal if a macOS sandbox blocks the watchdog's process inspection.

## Capabilities and conventions

`wavefunction_container="mps"` explicitly requests the optional Microsoft
`MPSContainer`/`MPSSite` API. It fails before solving if that capability is absent.
Small statevector examples work without it. `rdm` requires `RdmContainer`, which
is a separate capability; a QDK build may supply one without the other. The
`auto` policy is unchanged (small CI vector, otherwise RDM container).

The MPS bridge converts the doubly occupied local-state phase and uses lossless
QR factorizations to move the orthogonality center to site zero. It preserves
the final state and its spin-resolved RDMs. There is no determinant expansion,
extra optimization sweep or additional truncation. Export currently uses dense
site tensors with trivial symmetry labels: tensor storage is O(L m²), while
the integral and full 2-RDM payload remains O(L⁴). It does not demonstrate
symmetry-block compression in Microsoft's circuit constructor.

Reported calculator energy belongs to the returned final state, including its
core energy. Intermediate Lanczos energies and discarded weights remain
diagnostics. Adjacent iteration energies are not a completed-sweep convergence
test. The four-site spin-1 evolution example keeps all sites explicit in the
classic superblock, so it tests time-step convergence, not a bond-truncation
comparison. The six-site spin-1/2 example also compares retained-state budgets.

## Scope and provenance

Adapted from the advanced examples in the private `skilledwolf/alpsx` worktree
`98ed/ALPS`, based on commit `1adba871e`, including its local nitrogen/scaling
edits. Each adapted file identifies its source and SHA-256. See
`source-manifest.json` for the exact snapshots. The port uses public
`pyalps.models`, `pyalps.ed` and `pyalps.classic_dmrg` interfaces.

Source-only MPO compilation/restart, one-site polishing and S²-penalty solver
paths are not exposed by this branch's molecular calculator and are omitted.
The small-CAS spin check is consolidated into `nitrogen/run.py solve`; it
measures S² independently and does not claim singlet targeting. The two source
report generators, which embedded backend-specific results, are replaced by
one [render.py](render.py) for fresh nitrogen/scaling JSON. Historical timings,
figures and large binary results are not copied into this branch.

## Executed validation

See [the validation snapshot](validation/README.md) for fresh results, tested
dependency configurations and observed resource limits.
