# CT-HYB-SEGMENT

`ct-hyb-segment` is the segment-picture hybridization-expansion continuous-time
quantum Monte Carlo impurity solver described in

> H. Hafermann, P. Werner, E. Gull, *Efficient implementation of the
> continuous-time hybridization expansion quantum impurity solver*,
> Comput. Phys. Commun. **184**, 1280 (2013).

It was developed as [ALPSCore/CT-HYB-SEGMENT](https://github.com/ALPSCore/CT-HYB-SEGMENT)
on top of the ALPSCore libraries and is now maintained here, built against the
ALPS libraries. Its history is part of this repository. `ct-hyb-segment --citations`
lists the papers to cite.

It descends from the same code as the `hybridization` solver in
`../hybridization`, and the two are built side by side. They produce the same
results; see the conventions below.

## Building

The solver is built with the ALPS applications when Eigen 3.3 or newer is found
(it uses `ALPS::gf`). Options:

- `ALPS_ENABLE_MPI` — run on all ranks of `MPI_COMM_WORLD`.
- `ALPS_ENABLE_NFFT` (default `OFF`) — evaluate frequency measurements with
  [NFFT3](https://www-user.tu-chemnitz.de/~potts/nfft/) (set `NFFT3_DIR`); without
  it they are summed directly.

## Running

```
ct-hyb-segment parameters.ini [NAME=value ...]   # results in sim.h5
ct-hyb-segment input.h5 output.h5                # as the DMFT loop calls it
ct-hyb-segment --help                            # lists all parameters
```

- From the DMFT program, set `SOLVER = ct-hyb-segment` (with `SC_WRITE_DELTA`);
  parameters are translated as for `hybridization`.
- From Python, `pyalps.cthyb_segment.solve(parms)` takes the parameters below.
  `pyalps.cthyb_segment_compat.solve(parms)` takes parameters written for
  `pyalps.cthyb` and translates them.
- The hybridization and DMFT tutorials run on this solver when the environment
  variable `ALPS_CTHYB_SOLVER=ct-hyb-segment` is set.

A minimal parameter file for an impurity with two flavors:

```
BETA=10
U=2
MU=0                  # relative to half filling, see below
FLAVORS=2
N=500                 # imaginary-time points
NMATSUBARA=100
MAX_TIME=60
cthyb.DELTA="delta.dat"   # index and Delta(tau) per flavor, N+1 lines
cthyb.N_MEAS=10
cthyb.SWEEPS=1000000
cthyb.THERMALIZATION=10000
```

## Conventions that differ from `hybridization`

| | `ct-hyb-segment` | `hybridization` |
|---|---|---|
| solver options | prefixed, e.g. `cthyb.SWEEPS`, `cthyb.MEASURE_freq` | `SWEEPS`, `MEASURE_freq` |
| grid sizes, flavors | `N`, `NMATSUBARA`, `FLAVORS` | `N_TAU`, `N_MATSUBARA`, `N_ORBITALS` |
| a sweep | one update | `N_MEAS` updates |
| `MU` | relative to half filling: the solver adds `U/2` | absolute: half filling at `MU=U/2` |
| inter-orbital interaction | `Uprime`, default 0; Hund's `cthyb.J` | `U'`, default `U-2J`; `J` |
| `U_MATRIX` text file | dense `FLAVORS`×`FLAVORS` matrix | dense matrix |
| G(iω_n) measured | by default (`cthyb.MEASURE_freq=1`) | on request |

`MU_VECTOR` gives absolute orbital energies in both. With a retarded
interaction (`cthyb.RET_INT_K`), both solvers shift `U → U − 2K'(0)` and
`μ → μ + K'(0)`, as described in `Documentation/hybdoc.tex`.

`Documentation/hybdoc.tex` documents the algorithm, the input files and the
measurements; its parameter names are those of `hybridization`.
