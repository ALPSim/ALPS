# Frozen application results

These are complete task-result archives from the pre-modularization ALPS revision
`841fba675996a99b7adacc113c6d5852cb69d978` (version 3.0.0 development), generated
with the adjacent inputs. The existing native SDK build used GCC and the shared
Boost/HDF5 dependencies. `expected.json` records values read by that revision's
pyalps installation, independently of the current reader.

For each input, in a fresh directory with that SDK on PATH:

```sh
cp <application>.parm input
parameter2xml input
<executable> --Tmin 1 --Tmax 1 --write-xml input.in.xml
```

Executables are `spinmc`, `loop` (looper), `fulldiag` and `dmrg`.
Each `.h5.gz` losslessly compresses the corresponding `input.task1.out.h5`
with a zero gzip timestamp; no HDF5 datasets were rewritten. Tests unpack it
locally and exercise `pyalps.loadMeasurements` or `loadEigenstateMeasurements`.
They cover means/errors/counts, spectra and a local vector observable, including
multiple full-diagonalization sectors. They do not assert stochastic simulation
convergence: scheduling can change Monte Carlo sample counts on regeneration.

Keep these files frozen when replacing serializers or readers. A current-writer
round trip cannot replace this compatibility check. These files establish ALPS
result compatibility, not ALPSCore writer compatibility or restart coverage.
