# pyalps

Python applications and libraries for the Algorithms and Libraries for Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Matplotlib plotting helpers are included with `pyalps`. Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

Wheels bundle simulation applications for `pyalps.runApplication`, but do not add them to your shell's `PATH` or include auxiliary tools such as `parameter2xml`, `printgraph`, and `alps-xml`. Install the [C++ SDK and tools](../../CONTRIBUTING.md#build) for those command-line workflows. Tutorial files are available in the [source collection](../../tutorials/README.md), with an optional SDK installation component.

The bindings are built as a standalone `scikit-build-core` project using nanobind. A source build requires GIL-enabled CPython 3.10 or newer, CMake 3.27 or newer, Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed shared ALPS C++ SDK. The SDK's numeric version must match `ALPS_VERSION.txt`; CMake rejects a mismatch before compiling. Point `ALPS_DIR` at the SDK's `share/alps` package directory.

Before building the SDK below, follow the [CMake and Ninja setup](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system CMake is older than 3.27.

The `distribution` CMake preset builds the SDK exactly as the wheel CI does. From the repository root:

```sh
cmake --preset distribution
cmake --build --preset distribution

export ALPS_DIR="$PWD/_build/distribution/install/share/alps"
python -m pip install build
python -m build --wheel python/pyalps
```

The wheel is written to `python/pyalps/dist` and can be installed with `python -m pip install`. With ccache installed, configure with `cmake --preset distribution -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set `CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to speed up rebuilds.

`PYALPS_BUILD_SOLVERS=ON` is the default and builds the MaxEnt, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake configuration for a smaller core-only developer build. These modules link `ALPS::maxent`, `ALPS::cthyb`, and `ALPS::ctint` from an SDK built with `ALPS_BUILD_APPLICATIONS=ON`; they do not compile solver implementations.

All binding sources and build helpers live under this package directory. The source distribution includes shared release metadata from `cmake/ALPS_VERSION.txt` and `LICENSE.txt`, stored at its root as `ALPS_VERSION.txt` and `LICENSE.txt`. It can be built outside the checkout against an installed SDK, without the C++ source tree.

The [Ising extension example](examples/ising/README.md) demonstrates exporting a simulation through nanobind against the installed SDK and pyalps wheel.

`PYALPS_BUNDLE_APPLICATIONS=ON` independently controls inclusion of the ALPS command-line programs (`spinmc`, `dmrg`, `sparsediag`, `loop`, `qwl`, ...) in `pyalps/bin`. Required shared libraries are included even when this option is `OFF`. The `runApplication` helpers first honor programs already on `PATH`, then search `ALPS_BIN_PATH` and the bundled directory. A bindings-only installation can therefore use separately installed programs through either environment variable. XML resources come from the package unless `ALPS_XML_PATH` is set; installations contain no fallback paths to the machine that built the wheel.

## Editable development

Build and install the C++ SDK using [CMake](../../CONTRIBUTING.md#build), then use the Python environment of your choice for an editable installation. Keep the SDK and bindings on the same compiler, architecture and native dependency stack.

After installing the matching SDK, install build dependencies in your active Python environment and use pip's editable mode:

```sh
export ALPS_DIR="$PWD/_build/distribution/install/share/alps"
python -m pip install "scikit-build-core>=1.0" "nanobind==2.15.0" "cmake>=3.27" ninja \
  "patchelf>=0.14; sys_platform == 'linux'"
python -m pip install --no-build-isolation -e python/pyalps \
  --config-setting "build-dir=$PWD/_build/manual-python"
```

The example uses the SDK installed by the `distribution` preset above. For another SDK, change `ALPS_DIR` to its installed `share/alps` directory. The Linux `patchelf` dependency is required when bundling applications; `--no-build-isolation` makes installing build dependencies your responsibility. Pip installs NumPy, SciPy, and Matplotlib as runtime dependencies.

For a smaller SDK and core-only editable Python installation on Linux/macOS, use the `sdk` preset and disable both solver bindings and bundled programs. After installing the build dependencies above:

```sh
cmake --preset sdk
cmake --build --preset sdk --parallel 2
cmake --install _build/sdk
export ALPS_DIR="$PWD/_build/sdk/install/share/alps"
python -m pip install --no-build-isolation -e python/pyalps \
  --config-setting "build-dir=$PWD/_build/manual-python-core" \
  --config-setting cmake.define.PYALPS_BUILD_SOLVERS=OFF \
  --config-setting cmake.define.PYALPS_BUNDLE_APPLICATIONS=OFF
```

Python edits take effect in a new interpreter without reinstalling. Rebuild and reinstall the SDK after SDK C++ changes, then rerun the editable install after changing binding C++ sources, build configuration, or packaged runtime resources. If you switch SDKs, compilers, or dependency providers while reusing a manual binding build directory, add `--config-setting cmake.args=--fresh` once to clear cached discovery. Scikit-build-core owns the Python source mapping; CMake installs the native modules, libraries, XML, downstream headers, and CMake package. Resource lookup uses the installed runtime manifest, so it also works when these files live separately from the Python sources.

## Native runtime layout

The extensions share one nanobind library, named `pyalps_nanobind` to avoid filename collisions with other packages, and one shared ALPS runtime. On Unix, CMake installs relative runtime paths with the package's `lib` directory first and derives any additional dependency directories from resolved link targets. Linux and macOS wheels intended for redistribution must then be repaired with auditwheel or delocate, respectively, to bundle external dependencies and replace build-machine paths; the wheel CI performs this step. See [downstream native extensions](#downstream-native-extensions) for manifest finalization after repair.

## Python compatibility

pyalps supports GIL-enabled CPython 3.10 and newer. Wheels use the ABI of the interpreter that built them; they do not use the limited API or `abi3`. Free-threaded Python is unsupported.

Build dependencies pin nanobind 2.15.0. The installed `runtime.json` records both its version and internals ABI. The downstream CMake package checks that ABI before creating `pyalps::runtime`, so an incompatible nanobind installation produces a configuration error rather than an interpreter abort. Rebuild native consumers against the same SDK and nanobind release as the installed pyalps package.

## Versioning

pyalps does not carry a version of its own. The numeric version is read from `cmake/ALPS_VERSION.txt` — the same file `cmake/ALPSVersion.cmake` reads for `ALPS_VERSION_CORE` — so a release bump is one edit rather than two that can drift. `tests/pyalps/test_wheel_payload.py` fails if the installed version and that file disagree.

A prerelease label cannot live in that file: `project(VERSION ...)` rejects a non-numeric version, and neither `find_package()` matching nor the library SOVERSION has a notion of prerelease ordering. CMake takes it from the `ALPS_VERSION_PRERELEASE` cache variable; the Python build takes it from the environment variable of the same name, using the same vocabulary:

| `ALPS_VERSION_PRERELEASE` | version with `ALPS_VERSION.txt` = 2.3.4 |
|---|---|
| unset | `2.3.4` |
| `beta.1` | `2.3.4b1` |
| `alpha.2` | `2.3.4a2` |
| `rc.1` | `2.3.4rc1` |
| `dev.3` | `2.3.4.dev3` |

In GitHub release builds, the provider takes the prerelease label from `GITHUB_REF` (for example, `refs/tags/v3.0.0-beta.1`). It rejects a tag whose numeric version differs from `ALPS_VERSION.txt`, or whose label conflicts with an explicit `ALPS_VERSION_PRERELEASE`. Wheels and source distributions use the same provider. `python python/pyalps/_build_support/alps_version.py` prints the version a build would produce, from any working directory. An sdist preserves its recorded version when rebuilt without the original build environment.

Note the consequence: because the number is inherited, a Python-only API change cannot be signalled in the pyalps version alone — it takes a bump of `cmake/ALPS_VERSION.txt`, which moves the whole project.
## Downstream native extensions

The C++ SDK supplies `ALPS::alps` for standalone programs. The installed Python package separately supplies `pyalps::runtime` for extensions that share ALPS objects or HDF5 handles with pyalps. A matching C++ SDK is still required for headers and compile settings.

```cmake
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
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
nanobind_add_module(my_module NB_STATIC my_module.cpp)
target_link_libraries(my_module PRIVATE pyalps::runtime)
```

`pyalps.get_cmake_dir()` exposes the same directory to Python tools. The C++ SDK neither installs this package nor discovers Python. The former SDK function `alps_target_link_pyalps` has been removed.

The same target supplies `<pyalps/export_simulation.hpp>` for exporting a derived simulation through nanobind. This Python-owned header replaces the old SDK header `<alps/ngs/detail/export_sim_to_python.hpp>`; update that include when rebuilding a downstream extension.

Wheel installation writes `pyalps/runtime.json`. After auditwheel or delocate repair, regenerate it with `python python/pyalps/_build_support/runtime_manifest.py --wheel path/to/pyalps.whl`. Cibuildwheel runs this automatically. The manifest records the final relative library paths; pyalps exposes these as imported CMake targets. On macOS, wheel finalization sets linkable `@rpath` library IDs and refreshes their signatures, so downstream builds need no binary-patching commands.

CMake derives build-time search paths from the imported targets. If you install or redistribute your extension, set its `INSTALL_RPATH` for the destination layout using normal CMake installation rules. The target does not hard-code the build environment's Python installation into installed extensions. For an extension installed for the same environment, CMake's `INSTALL_RPATH_USE_LINK_PATH` target property can retain the runtime search paths.
