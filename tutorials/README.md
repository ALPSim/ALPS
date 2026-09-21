# Learn ALPS

These tutorials are the starting point for learning ALPS: prepare a simulation, run it, inspect its results, and then explore a numerical method or build your own code. The [library examples](examples/README.md) provide smaller demonstrations of individual APIs.

## Getting started

1. Follow the [installation instructions](../README.md#installation) to install ALPS and its Python package, `pyalps`.
2. Work through [the basics](intro-01-basics/). Start with `tutorial-prepareinput.py`, `tutorial-runsimulation.py`, and `tutorial-evaluate.py`; `tutorial-full.py` combines the workflow.
3. Continue with [autocorrelations](mc-01-autocorrelations/) and [equilibration and convergence](mc-01b-equilibration-and-convergence/) before interpreting Monte Carlo measurements.

Run scripts from their own tutorial directory using the Python environment in which `pyalps` is installed. Starting from this directory:

```sh
cd intro-01-basics
python tutorial-full.py
```

For an interactive presentation, browse the [English notebooks](notebook/en/) or [Japanese notebooks](notebook/ja/).

## Running simulations

Choose a method, then follow its numbered tutorial sequence.

| Method or task | Suggested starting point | Continue with |
| --- | --- | --- |
| Classical Monte Carlo | [Autocorrelations](mc-01-autocorrelations/) | [Susceptibilities](mc-02-susceptibilities/), [magnetization](mc-03-magnetization/), [measurements](mc-04-measurements/) |
| Quantum Monte Carlo | [Bosons](mc-05-bosons/) | [Quantum Wang–Landau](mc-06-qwl/), [quantum phase transitions](mc-08-quantum-phase-transition/) |
| Exact diagonalization | [Sparse diagonalization](ed-01-sparsediag/) | [Gaps](ed-02-gaps/), [spectra](ed-03-1dspectra/), [full diagonalization](ed-06-fulldiag/) |
| DMRG | [Ground-state energies](dmrg-03-ground-state-energies/) | [Gaps](dmrg-04-gaps/), [local observables](dmrg-05-local-observables/), [correlations](dmrg-06-correlations/) |
| DMFT | [Hybridization expansion](dmft-02-hybridization/) | [Interaction expansion](dmft-03-interaction/), [Mott transition](dmft-04-mott/), [orbital-selective transitions](dmft-05-osmt/) |
| Python impurity-solver interface | [Hybridization solver](hybridization-01-python/) | [Kondo model](hybridization-02-kondo/), [retarded interactions](hybridization-03-retarded-interaction/), [spin freezing](hybridization-04-spinfreezing/) |

## Solver references

- [Looper](looper/README.md): path-integral and SSE algorithms, parameters, annealing, and measurements.
- [Hybridization expansion](hybridization/README.md): running CT-HYB and its detailed scientific manual.

## Developing with ALPS

- Start with the [Python simulation skeleton](code-01-python/) or [C++ simulation example](code-02-c++/).
- Follow the `alpsize-*` sequence from [the CMake introduction](alpsize-01-cmake/) through parameters, measurements, lattices, and scheduling.
- Explore the [accumulator example](ngs/1_accumulator_only/) and the [native Python simulation](ngs/6_python_native/) for simulation interfaces.
- Export a C++ simulation to Python with the [Ising extension example](../python/pyalps/examples/ising/README.md), built against the installed SDK and pyalps package.
- Use the [library examples](examples/README.md) when you need a focused example of statistical analysis, HDF5, model construction, numerical methods, or Fortran integration.

## Using an installed tutorial collection

Source builds provide this collection as an optional installation component:

```sh
cmake --install build --component tutorials
```

It is installed under `share/alps/tutorials`. The library examples can also be built against an installed SDK; see their [build instructions](examples/README.md).
