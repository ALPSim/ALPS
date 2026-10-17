# pyalps

Python applications and libraries for the [ALPS project](https://alps.comp-phys.org/). Install from PyPI in a virtual environment:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install pyalps
```

Use a separate environment to avoid overwriting commands from a source-installed ALPS SDK, or use the bindings-only installation described below.

Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` interface. Bundled applications are serial even with this extra installed.

## Command-line tools

Wheels include commands for the bundled applications, including `spinmc`, `loop`, `worm`, `dmrg`, and `sparsediag`, plus the tutorial tools `parameter2xml`, `printgraph`, `convert2xml`, `convert2text`, `plot2text`, `plot2gp`, `plot2xmgr`, `snap2vtk`, and `maxent`. For example:

```sh
parameter2xml simulation.in
spinmc --write-xml simulation.in.in.xml
convert2text simulation.in.task1.out.xml > results.txt
```

Download inputs from the [ALPS tutorials](https://alps.comp-phys.org/tutorials/); they are not included in wheels. Set `ALPS_XML_PATH` only to override the bundled XML resources.

For older tutorials, use `python` instead of the removed `alpspython` wrapper and `pyalps.plot` instead of `plot2mpl` or `extractmpl`.

## Building from source

Requires GIL-enabled CPython 3.10+, CMake 3.27+, Ninja, a C++17 compiler, external Boost 1.76+, LP64 BLAS/LAPACK, and HDF5 1.10.5+. Reuse an installed ALPS C++ SDK, or build one with the `distribution` preset. From the repository root, build the SDK and Python wheel with:

```sh
cmake --preset distribution
cmake --build --preset distribution

python -m pip install build
ALPS_DIR="$PWD/_build/distribution/install/share/alps" \
  python -m build --wheel python/pyalps
```

Install the wheel from `python/pyalps/dist` with `python -m pip install`. It bundles the SDK's applications and libraries; shell launchers and Python helpers use those executables regardless of `PATH` or `ALPS_BIN_PATH`.

For bindings only, set `ALPS_DIR` to the existing SDK's `share/alps` directory and set the environment variable `PYALPS_BUNDLE_APPLICATIONS=OFF`:

```sh
ALPS_DIR="/path/to/sdk/share/alps" PYALPS_BUNDLE_APPLICATIONS=OFF \
  python -m pip install ./python/pyalps
```

This installs no command launchers and can safely share the SDK's prefix. Keep the SDK installed and add its `bin` to `PATH` for shell use. Python helpers use `ALPS_BIN_PATH` or the SDK recorded in the installed runtime manifest; they do not search `PATH`.

In either mode, pass a full executable path to a Python helper to select a different application, including a source-built MPI application.

## Compatibility

- Python 3.10 and 3.11 have separate `cp310` and `cp311` wheels; Python 3.12+ shares a `cp312-abi3` wheel. Free-threaded Python is unsupported.
- Rebuild downstream C++ extensions against the same SDK revision, nanobind internals ABI, and C++ standard library as the wheel. Nanobind is pinned to 2.15.0; its version and internals ABI are recorded in the installed `pyalps/runtime.json`.
- Legacy HDF5 signed-byte datasets without type metadata are read as booleans. Use a typed reader such as h5py when such a dataset contains integers.
- `pyalps.mpi` uses mpi4py's protocol, which is incompatible with Boost.MPI serialization. Communicators cannot be used as dictionary keys.

The package version comes from `cmake/ALPS_VERSION.txt`. For local prereleases, set `ALPS_VERSION_PRERELEASE` (for example, `beta.1`); release builds derive it from the Git tag.

## Editable development

After installing the matching SDK and build dependencies, use `python -m pip install --no-build-isolation -e python/pyalps`. Python edits take effect in a new interpreter; rebuild after native changes. Installed runtime resources are located through `runtime.json`, separately from editable Python sources. See [contributor setup](../../CONTRIBUTING.md#build).

`PYALPS_BUILD_SOLVERS=OFF` disables MaxEnt, CT-HYB and CT-INT bindings independently of bundling executables. Pass it as `--config-setting cmake.define.PYALPS_BUILD_SOLVERS=OFF`; use the `PYALPS_BUNDLE_APPLICATIONS=OFF` environment setting when building against an SDK without applications. Solver bindings otherwise link the installed SDK libraries without recompiling their implementations.

## Models, exact diagonalization and classic DMRG

`pyalps.models`, `pyalps.ed` and `pyalps.classic_dmrg` bind the SDK components [`ALPS::models`](../../src/alps/models/README.md), [`ALPS::ed`](../../src/alps/ed/README.md) and [`ALPS::classic_dmrg`](../../src/alps/classic_dmrg/README.md): lattices and models built in Python or imported from the XML catalog, compiled to numeric local terms, and solved by exact diagonalization or DMRG without the XML applications.

```python
from pyalps import models, ed, classic_dmrg as dmrg
lat = models.lattice

model = models.heisenberg(lat.chain(8))                    # periodic ring, J = 1
problem = ed.prepare_product(model, sector={"Sz": 0})      # one charge sector, 70 states
energies = problem.eigenstates(count=2).eigenvalues        # dense or Lanczos by block size
spectrum = ed.model_spectrum(model, all_momenta=True)      # every (Sz, k) block

spin = models.spin()
custom = models.Hamiltonian(lat.square(4, 4), spin, conserved_quantum_numbers=["Sz"])
custom.add_bond(1.0, lambda i, j: spin.Sz(i) * spin.Sz(j))     # spin.Sz(i) is models.located(spin, "Sz", i)
ground = ed.ground_state(model, sector={"Sz": 0},
                         observables={"SzSz": spin.Sz(0) * spin.Sz(1), "Sz": "Sz"})
dense = model.to_dense_matrix(sector={"Sz": 0})             # or to_sparse_matrix() for scipy CSR

pauli = models.export_pauli(model)                         # Pauli strings for qubit simulators
catalog = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 8, "J": 1,
                               "CONSERVED_QUANTUMNUMBERS": "Sz", "Sz_total": 0})
levels = ed.prepare_product(catalog).eigenstates(count=2).eigenvalues  # its Sz_total sector

chain = models.hubbard(lat.chain(16, bc="open"), t=1.0, U=4.0)
result = dmrg.solve(chain, sweeps=dmrg.ramped_sweeps(4, 200), sector={"N": 16, "Sz": 0},
                    measurements={"n": dmrg.site_terms(chain, "n")})
```

The `ed` entry points (`prepare_product`, `prepare_projected`, `prepare_finite`, `model_spectrum`, `finite_spectrum`, `ground_state`) and `dmrg.solve` take a kit `Hamiltonian` or an XML `CatalogModel` alike; a catalog runs in the sector of its `*_total` parameters unless `sector=` says otherwise, and resolves momenta only with explicit `translations=`. Observables (`ed.ground_state(observables=...)`, product results' `expectation_values`, `dmrg.solve(measurements=..., expectation=...)`) may be compiled operators, located terms or operator names; `models.compile_observable` does the compiling and rejects terms built on another local basis, such as kit spins on a catalog model, whose Sz states run the other way. Site operators answer to either the kit's or the XML catalog's names (`Sp`/`Splus`, `c_dag_up`/`cdag_up`, `n_dn`/`n_down`; `models.catalog_operator_names`), and `TermSum.dag()` names conjugates as the site does (`Sp` → `Sm`).

Coefficients may be symbols: `J, h = models.symbols("J", "h")` and `J * spin.Sz(0) * spin.Sz(1)` is a `SymbolicTermSum`; a `models.ParametricHamiltonian` takes the `Hamiltonian` arguments and symbolic `add_bond`/`add_site`/`add` coefficients, and `bind(J=1.0, h=0.1)` returns the numeric `Hamiltonian` the solvers take (a partial bind stays parametric). `models.validate(model, target=None | "ed" | "dmrg")` returns a report of structural and solver-specific diagnostics (conservation, Hermiticity, unbound symbols, sector size, lattice shape), and `models.has_sign_problem(model)` runs the legacy algebraic QMC sign check on a kit or catalog model.

`warmup_configuration=[...]` (one local-state index per site, in each site's `PhysicalBasis` order) sets the charges the DMRG warm-up targets on each partial chain to that product configuration's, instead of uniform shares of the sector. Ground states with unevenly distributed charges need it, such as a domain wall from fields alone or random fields: from a uniform warm-up their sweeps keep misplaced charges and report a wrong energy without warning ([classic DMRG README](../../src/alps/classic_dmrg/README.md#unevenly-distributed-ground-states)). `dmrg.evolve` passes it to the preparation of the initial state.

`keep_mps=True` also returns every targeted state as a `dmrg.MatrixProductState` in `result.mps`: the engine's final block rotations and closing wave function, so exactly the state the energies and measurements belong to. Its `tensors` are NumPy arrays `[left bond, local state, right bond]` in each site's `PhysicalBasis` order, mixed canonical around `center`; `to_dense()`, `amplitude(states)`, `overlap(other)` and `norm()` contract them, and on small chains `to_dense()` is the `pyalps.ed` eigenvector up to a phase. The XML `dmrg` application writes the same state with `SAVE_MPS=1`, read by `pyalps.loadMPS`.

`dmrg.evolve(initial, final, time_step=..., steps=..., sector=...)` evolves the DMRG ground state of `initial` under `final` by adaptive time-dependent DMRG: second-order Trotter gates on the bonds, applied while sweeping. Real time runs in complex arithmetic and `imaginary=True` applies `exp(-tau H)`; `final` may couple only neighbouring sites in chain order, and `drives=[(operator, f)]` add `f(t) * operator` (sampled at the step midpoints). `method="targeting"` uses time-step targeting with a Krylov step instead: `final` may couple any sites, there is no time-step error and there are no drives; `half_sweeps` (2) sets the half sweeps per step. Every measured step (`measure_steps`) reports its time, `<H_final>` and the measurements, with `keep_mps=True` its state; `dmrg.overlaps` and `dmrg.loschmidt_echo` contract the kept states. Either Hamiltonian may be an operator (a `TermSum` or `CompiledOperator`) on the other's basis. [tutorials/time-evolution](../../tutorials/time-evolution/README.md) has quench, domain-wall and Loschmidt-echo examples.

`dmrg.solve_molecular(h1, eri, n_alpha, n_beta, core_energy=..., sweeps=...)` runs the same engine on restricted real molecular integrals (chemist notation) and returns the energy with the one- and two-particle reduced density matrices in the PySCF convention, spin-summed and by spin block. Everything is measured in the final sweep: `observables={name: TermSum}` adds operators of up to four orbitals (`models.located(models.spinful_fermion(), op, orbital)` products), and `ci_vector=True` reconstructs the sector's CI vector (`dmrg.ci_occupations(result)` gives its alpha and beta occupations); `keep_mps=True` returns the final state as an MPS over the orbitals. `dmrg.local_expectations`, `dmrg.correlation_matrix` and `dmrg.spin_spin_correlations` derive densities, number, spin and hopping correlations from the RDMs. The solver suits small and medium active spaces, up to about 16 orbitals (see the [component README](../../src/alps/classic_dmrg/README.md#molecular-hamiltonians-alpsdmrgmolecularhpp)). [`qdk-chemistry-alps`](../qdk-chemistry-alps/README.md) registers exact diagonalization and this solver as qdk-chemistry calculators.

Matrices and spectra come back as NumPy arrays (eigenvectors column-major; sparse matrices convert with `to_scipy()`), and the compile and solve calls release the GIL. The module docstrings and `help()` describe the Python surface; the C++ READMEs document the conventions (fermion sign conventions, momentum labels, kept-state schedules). Bond and site specifications are Python callables, called once per bond or site when added.

`PYALPS_BUILD_MODELS=OFF` (`--config-setting cmake.define.PYALPS_BUILD_MODELS=OFF`) skips these bindings. `pyalps.models` and `pyalps.ed` need only the always-built `models` and `ed` SDK components; `pyalps.classic_dmrg` is present when the SDK was built with `ALPS_ENABLE_CLASSIC_DMRG` (the default, which needs Eigen) and raises `ImportError` otherwise.

## Downstream native extensions

The C++ SDK supplies `ALPS::alps` for standalone programs. The installed Python package separately supplies `pyalps::runtime` for extensions that share ALPS objects or HDF5 handles with pyalps. A matching C++ SDK is still required for headers and compile settings.

```cmake
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE nanobind_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind 2.15.0 EXACT CONFIG REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pyalps --cmake-dir
  OUTPUT_VARIABLE pyalps_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(pyalps CONFIG REQUIRED)
nanobind_add_module(my_module STABLE_ABI NB_STATIC my_module.cpp)
target_link_libraries(my_module PRIVATE pyalps::runtime)
```

`pyalps.get_cmake_dir()` exposes the same directory to Python tools. The C++ SDK neither installs this package nor discovers Python. The former SDK function `alps_target_link_pyalps` remains available, deprecated, from `include(${ALPS_PYTHON_USE_FILE})`; it requires `find_package(nanobind)` first and links `pyalps::runtime`.

The same target supplies `<pyalps/export_simulation.hpp>` for exporting a derived simulation through nanobind. The historical SDK include `<alps/ngs/detail/export_sim_to_python.hpp>` forwards to this Python-owned header and remains available.

Wheel installation writes `pyalps/runtime.json`. After auditwheel or delocate repair, regenerate it with `python python/pyalps/_build_support/runtime_manifest.py --wheel path/to/pyalps.whl`. Cibuildwheel runs this automatically. The manifest records the final relative library paths; pyalps exposes these as imported CMake targets. On macOS, wheel finalization sets linkable `@rpath` library IDs and refreshes their signatures, so downstream builds need no binary-patching commands.

CMake derives build-time search paths from the imported targets. If you install or redistribute your extension, set its `INSTALL_RPATH` for the destination layout using normal CMake installation rules. The target does not hard-code the build environment's Python installation into installed extensions. For an extension installed for the same environment, CMake's `INSTALL_RPATH_USE_LINK_PATH` target property can retain the runtime search paths.

The historical Boost.Python comparison harness has been retired. See [the coverage assessment](../../tests/python-migration.md) for its replacement regressions and the remaining limits of historical compatibility.
