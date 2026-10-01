# Learn ALPS

These tutorials are the starting point for learning ALPS: prepare a simulation, run it, inspect its results, and then explore a numerical method or build your own code. The [library examples](00-examples/README.md) provide smaller demonstrations of individual APIs.

## Getting started

1. Follow the [installation instructions](../README.md#installation) to install ALPS and its Python package, `pyalps`.
2. Work through [the basics](01-intro/). Start with `tutorial-prepareinput.py`, `tutorial-runsimulation.py`, and `tutorial-evaluate.py`; `tutorial-full.py` combines the workflow.
3. Continue with [autocorrelations](03-mc/01-autocorrelations/) and [equilibration and convergence](03-mc/01b-equilibration-and-convergence/) before interpreting Monte Carlo measurements.

Run scripts from their own tutorial directory using the Python environment in which `pyalps` is installed. Starting from this directory:

```sh
cd 01-intro
python tutorial-full.py
```

For an interactive presentation, browse the [English notebooks](11-notebook/en/) or [Japanese notebooks](11-notebook/ja/).

## Running simulations

Choose a method, then follow its numbered tutorial sequence.

The top-level numbers suggest a browsing order; each topic keeps its existing lesson numbers. Choose the methods relevant to your work rather than completing every group.

| Directory | Purpose |
| --- | --- |
| `00-examples/` | Independent API examples to consult as needed |
| `01-intro/` | Prepare, run and evaluate a first simulation |
| `02-ed/` | Exact diagonalization and small-system calculations |
| `03-mc/` | Classical and quantum Monte Carlo |
| `04-dmrg/` | Density-matrix renormalization group |
| `05-dmft/` | Dynamical mean-field theory |
| `06-hybridization/` | Impurity-solver tutorials and reference |
| `07-looper/` | Looper solver reference |
| `08-alpsize/` | Progressively integrate a simulation with ALPS |
| `09-code/` | Simulation skeletons and implementation examples |
| `10-ngs/` | Accumulator and simulation-interface examples |
| `11-notebook/` | Interactive versions of lessons, usable alongside the other topics |

`00-examples` is a reference collection, not a prerequisite. For writing simulations, follow `08-alpsize` before using the templates in `09-code` and the interfaces in `10-ngs`. Notebooks provide an alternative presentation rather than a final advanced lesson.

| Method or task | Suggested starting point | Continue with |
| --- | --- | --- |
| Classical Monte Carlo | [Autocorrelations](03-mc/01-autocorrelations/) | [Susceptibilities](03-mc/02-susceptibilities/), [magnetization](03-mc/03-magnetization/), [measurements](03-mc/04-measurements/) |
| Quantum Monte Carlo | [Bosons](03-mc/05-bosons/) | [Quantum Wang–Landau](03-mc/06-qwl/), [quantum phase transitions](03-mc/08-quantum-phase-transition/) |
| Exact diagonalization | [Sparse diagonalization](02-ed/01-sparsediag/) | [Gaps](02-ed/02-gaps/), [spectra](02-ed/03-1dspectra/), [full diagonalization](02-ed/06-fulldiag/) |
| DMRG | [Ground-state energies](04-dmrg/03-ground-state-energies/) | [Gaps](04-dmrg/04-gaps/), [local observables](04-dmrg/05-local-observables/), [correlations](04-dmrg/06-correlations/) |
| DMFT | [Hybridization expansion](05-dmft/02-hybridization/) | [Interaction expansion](05-dmft/03-interaction/), [Mott transition](05-dmft/04-mott/), [orbital-selective transitions](05-dmft/05-osmt/) |
| Python impurity-solver interface | [Hybridization solver](06-hybridization/01-python/) | [Kondo model](06-hybridization/02-kondo/), [retarded interactions](06-hybridization/03-retarded-interaction/), [spin freezing](06-hybridization/04-spinfreezing/) |

## Solver references

- [Looper](07-looper/README.md): path-integral and SSE algorithms, parameters, annealing, and measurements.
- [Hybridization expansion](06-hybridization/README.md): running CT-HYB and its detailed scientific manual.

## Developing with ALPS

- Follow the [`08-alpsize/` sequence](08-alpsize/) from [the CMake introduction](08-alpsize/01-cmake/) through parameters, measurements, lattices, and scheduling.
- Adapt the [Python simulation skeleton](09-code/01-python/) or [C++ simulation example](09-code/02-c++/).
- Explore the [accumulator example](10-ngs/1_accumulator_only/) and the [native Python simulation](10-ngs/6_python_native/) for simulation interfaces.
- Export a C++ simulation to Python with the [Ising extension example](../python/pyalps/examples/ising/README.md), built against the installed SDK and pyalps package.
- Use the [library examples](00-examples/README.md) when you need a focused example of statistical analysis, HDF5, model construction, numerical methods, or Fortran integration.

## Using an installed tutorial collection

Source builds provide this collection as an optional installation component:

```sh
cmake --install build --component tutorials
```

It is installed under `share/alps/tutorials`. The library examples can also be built against an installed SDK; see their [build instructions](00-examples/README.md).
