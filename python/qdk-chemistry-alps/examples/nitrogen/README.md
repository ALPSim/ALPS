# Nitrogen bond breaking and QDK state preparation

`run.py prepare` constructs N₂/cc-pVDZ, uses AVAS to select N 2p orbitals,
optimizes a singlet CAS(6e,6o) with CASSCF and computes an independent PySCF FCI
reference. Four inactive doubly occupied orbitals enter the effective one-body
integrals and scalar core energy. The exported problem has 12 qubits.

From the parent `examples` directory, for each distance (1.1, 1.5, 2.0, 2.5 Å):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py prepare --distance 2.5 --folder results/n2-r2.5
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py solve --folder results/n2-r2.5 --bonds 8 16 32 64 --sweeps 12
python render.py results/n2-r2.5/solve.json
```

`solve` uses `alps_classic_dmrg` and exports the final MPS. It checks signed
amplitudes against an independent fermionic permutation, energy against its
RDM contraction, particle normalization, S², right-canonical isometries, and
HDF5 round trips. MACIS and PySCF supply independent energy references. Dense
vectors are limited to this six-orbital validation problem; they are not part
of the production MPS export. At stretched bonds, inspect fidelity and S² in
addition to energy: fixed Nα=Nβ is not a singlet constraint.

## Optional Microsoft PR workflow

The tested dependency snapshots are QDK Chemistry MPS-container PR
[#751](https://github.com/microsoft/qdk-chemistry/pull/751),
`079029fad1b706a0d59aa3f1ae603aa26e2ba913`, and state-preparation PR
[#540](https://github.com/microsoft/qdk-chemistry/pull/540),
`fb286eb2d9dde51d960f6953cb726cb45d52c12b` (includes that container snapshot).
These are development APIs; a released QDK version number alone does not
guarantee them. Use a compatible combined build for the following stages.

```bash
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py circuit --folder results/n2-r2.5 --bond 16 --bits 8
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py circuit --folder results/n2-r2.5 --bond 16 --bits 10
python scaling/guard.py --max-seconds 300 -- python nitrogen/run.py qpe_resources --folder results/n2-r2.5 --bond 16 --bits 8
```

`circuit` executes the preparation in Q#, keeps the virtual-bond register until
inspection and includes every bond outcome in the reduced system state. It
reports finite-precision leakage, fidelity and energy without postselection.
It saves synthesis parameters and uses the composable operation because the
pinned PR's standalone QIR wrapper resets the system register. Each resource
estimate uses a fresh Q# context. `--method block_sparse` also accepts the
trivial-sector tensors, but cannot exploit charge labels discarded by export.

`qpe_resources` defaults to a small 8-phase-bit, time-0.5, second-order,
4-division Trotter resource demonstration. It does not execute QPE. Its
49.1 mEh phase-bin width is deliberately coarse. The source scenario is
selected with `--phase-bits 14 --trotter-divisions 32` (0.767 mEh bins), but
its counter exceeded even a 16 GB physical-footprint guard in local validation.
Phase-bin width is not a total-error bound; Trotter, synthesis and physical-
hardware costs require separate analysis.

The 8-bit dense preparation completed locally; a 10-bit simulation exceeded
420 seconds. Use a longer time budget deliberately if running that precision.

Optional `run.py block2 --folder ...` uses a separately installed `block2`
Python wheel as a small-CAS reference. Execute it in its own process to avoid
mixing OpenMP runtimes. It is not required by the ALPS/QDK port.

## CAS(10e,26o) capacity experiment

```bash
python scaling/guard.py --max-seconds 300 -- python nitrogen/capacity.py prepare --distance 2.0 --folder results/n2-large
python scaling/guard.py --max-seconds 300 -- python nitrogen/capacity.py alps --folder results/n2-large --bond 32 --sweeps 4
```

Preparation freezes four core electrons and orders the 26 RHF orbitals with
an exchange-weighted Fiedler vector. It never calls a large CASCI solver.
There are 4,327,008,400 determinants; even one float64 CI vector needs about
34.6 GB. The example instead retains an MPS and full RDMs. The four-sweep
command is a capacity pilot, not a converged singlet calculation. Check S²,
repeat at larger retained-state budgets and more sweeps before drawing chemical
conclusions. There is no restart or singlet penalty option in this port.
