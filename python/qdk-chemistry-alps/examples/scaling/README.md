# Classic DMRG scaling examples

`run.py` solves open Hubbard chains or full-Coulomb H chains in STO-3G with
symmetrically orthogonalized atomic orbitals ordered along the chain. It
returns a QDK MPS (default) or an RDM container without constructing CI
coefficients. References are restricted to at most 50,000 determinants.

From the parent `examples` directory:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python scaling/guard.py --max-seconds 300 -- python scaling/run.py --case hydrogen --sites 8 --bonds 32 64 --sweeps 12 --output results/hydrogen8.json
python scaling/guard.py --max-seconds 300 -- python scaling/run.py --case hydrogen --sites 20 --bonds 32 --sweeps 4 --output results/hydrogen20.json
python scaling/guard.py --max-seconds 300 -- python scaling/run.py --case hubbard --sites 32 --bonds 32 --sweeps 4 --output results/hubbard32.json
python render.py results/hydrogen8.json
```

Run cases sequentially. The watchdog samples child-tree RSS and, on macOS,
physical footprint, stops above 12 GB or below 8 GB available system RAM by
default, and kills the process group on timeout/interruption. It is a sampled
guard, not an OS-enforced allocation limit. Process inspection may require
running outside an application sandbox.

JSON reports integral preparation, solver time including RDMs, cumulative
process peak RSS, dense MPS bytes, final-state energy/RDM/trace checks, natural
occupations and iteration diagnostics. RSS includes imports and any earlier
reference solve in the same process; use separate invocations for per-setting
memory measurements. These examples do not measure peak scratch usage or
establish asymptotic runtime scaling from a few cases.

The large commands are capacity checks. Neither a small final discarded weight
nor adjacent iteration energies certify convergence. Compare longer runs and
larger kept-state budgets. Integral screening has an explicit conservative
energy bound. The combinatorial CI-memory estimate is an analytical lower
bound, not an observed allocation failure. Full integrals and full 2-RDMs
remain O(L⁴), even though the MPS avoids the determinant expansion.
