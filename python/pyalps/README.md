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

Requires GIL-enabled CPython 3.11+, CMake 3.27+, Ninja, a C++17 compiler, external Boost 1.76+, LP64 BLAS/LAPACK, and HDF5 1.10.5+. Reuse an installed ALPS C++ SDK, or build one with the `distribution` preset. From the repository root, build the SDK and Python wheel with:

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

- Python 3.11 has a separate `cp311` wheel; Python 3.12+ shares a `cp312-abi3` wheel. Free-threaded Python is unsupported.
- Rebuild downstream C++ extensions against the same SDK revision, nanobind internals ABI, and C++ standard library as the wheel. Nanobind is pinned to 2.15.0; its version and internals ABI are recorded in the installed `pyalps/runtime.json`.
- Legacy HDF5 signed-byte datasets without type metadata are read as booleans. Use a typed reader such as h5py when such a dataset contains integers.
- `pyalps.mpi` uses mpi4py's protocol, which is incompatible with Boost.MPI serialization. Communicators cannot be used as dictionary keys.

## Versioning

The package version comes from `cmake/ALPS_VERSION.txt`. For local prereleases, set `ALPS_VERSION_PRERELEASE` (for example, `beta.1`); release builds derive it from the Git tag. The unreleased SDK changes follow the published 3.0 release; maintainers must select and update the next release version before tagging.

## Editable development

After installing the matching SDK and build dependencies, use `python -m pip install --no-build-isolation -e python/pyalps`. Python edits take effect in a new interpreter; rebuild after native changes. Installed runtime resources are located through `runtime.json`, separately from editable Python sources. See [contributor setup](../../CONTRIBUTING.md#build).

`PYALPS_BUILD_SOLVERS=OFF` disables MaxEnt, CT-HYB and CT-INT bindings independently of bundling executables. Pass it as `--config-setting cmake.define.PYALPS_BUILD_SOLVERS=OFF`; use the `PYALPS_BUNDLE_APPLICATIONS=OFF` environment setting when building against an SDK without applications. Solver bindings otherwise link the installed SDK libraries without recompiling their implementations.

## Downstream native extensions

The native component and solver interfaces are changes for a future release
after 3.0; see the [SDK boundaries and compatibility guidance](../../CONTRIBUTING.md#consuming-the-c-sdk).

The C++ SDK supplies `ALPS::alps` for standalone programs. The installed Python package separately supplies `pyalps::runtime` for extensions that share ALPS objects or HDF5 handles with pyalps. A matching C++ SDK is still required for headers and compile settings.

```cmake
find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module
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

`pyalps.get_cmake_dir()` exposes the same directory to Python tools. The C++ SDK neither installs this package nor discovers Python. The former SDK function `alps_target_link_pyalps` has been removed.

The same target supplies `<pyalps/export_simulation.hpp>` for exporting a derived simulation through nanobind. This Python-owned header replaces the old SDK header `<alps/ngs/detail/export_sim_to_python.hpp>`; update that include when rebuilding a downstream extension.

## Native runtime layout

Wheel installation writes `pyalps/runtime.json`. After auditwheel or delocate repair, regenerate it with `python python/pyalps/_build_support/runtime_manifest.py --wheel path/to/pyalps.whl`. Cibuildwheel runs this automatically. The manifest records the final relative library paths; pyalps exposes these as imported CMake targets. On macOS, wheel finalization sets linkable `@rpath` library IDs and refreshes their signatures, so downstream builds need no binary-patching commands.

CMake derives build-time search paths from the imported targets. If you install or redistribute your extension, set its `INSTALL_RPATH` for the destination layout using normal CMake installation rules. The target does not hard-code the build environment's Python installation into installed extensions. For an extension installed for the same environment, CMake's `INSTALL_RPATH_USE_LINK_PATH` target property can retain the runtime search paths.
