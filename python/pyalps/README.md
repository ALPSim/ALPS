# pyalps

Python applications and libraries for the Algorithms and Libraries for Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Install `pyalps[plot]` to use the Matplotlib plotting helpers. Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

The bindings are built as a standalone `scikit-build-core` project using nanobind. A source build requires GIL-enabled CPython 3.12 or newer, CMake 4.3 or newer, Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed shared ALPS C++ SDK. The SDK's numeric version must match `ALPS_VERSION.txt`; CMake rejects a mismatch before compiling. Point `ALPS_DIR` at the SDK's `share/alps` package directory.

Before building the SDK below, follow the [CMake and Ninja setup](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system CMake is older than 4.3.

The `distribution` CMake preset builds the SDK exactly as the wheel CI does. From the repository root:

```sh
cmake --preset distribution
cmake --build --preset distribution

ALPS_DIR="$PWD/_build/distribution/install/share/alps" \
  python -m build --wheel python/pyalps
```

The wheel is written to `python/pyalps/dist` and can be installed with `python -m pip install`. With ccache installed, configure with `cmake --preset distribution -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set `CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to speed up rebuilds.

`PYALPS_BUILD_SOLVERS=ON` is the default and builds the MaxEnt, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake configuration for a smaller core-only developer build. These modules link `ALPS::maxent`, `ALPS::cthyb`, and `ALPS::ctint` from an SDK built with `ALPS_BUILD_APPLICATIONS=ON`; they do not compile solver implementations.

All binding sources and build helpers live under this package directory. The source distribution includes only two files from the repository root: `ALPS_VERSION.txt` and `LICENSE.txt`, shared release metadata. It can be built outside the checkout against an installed SDK, without the C++ source tree.

The [Ising extension example](examples/ising/README.md) demonstrates exporting a simulation through nanobind against the installed SDK and pyalps wheel.

`PYALPS_BUNDLE_APPLICATIONS=ON` independently controls inclusion of the ALPS command-line programs (`spinmc`, `dmrg`, `sparsediag`, `loop`, `qwl`, ...) in `pyalps/bin`. Required shared libraries are included even when this option is `OFF`. The `runApplication` helpers first honor programs already on `PATH`, then search `ALPS_BIN_PATH` and the bundled directory. A bindings-only installation can therefore use separately installed programs through either environment variable. XML resources come from the package unless `ALPS_XML_PATH` is set; installations contain no fallback paths to the machine that built the wheel.

On Windows, first build and install the `windows-x64` or `windows-arm64` preset described in [CONTRIBUTING.md](../../CONTRIBUTING.md#native-windows-msvc). Use CPython and a Visual Studio developer shell targeting the same architecture as the SDK, then:

```powershell
$arch = 'arm64' # Use 'x64' for the x64 SDK and Python.
$env:ALPS_DIR = "$pwd/_build/windows-$arch/install/share/alps"
python -m pip install build
python -m build --wheel python/pyalps `
  --config-setting "cmake.define.CMAKE_TOOLCHAIN_FILE=$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  --config-setting "cmake.define.VCPKG_INSTALLED_DIR=$pwd/_build/windows-$arch/vcpkg_installed" `
  --config-setting cmake.define.VCPKG_MANIFEST_MODE=OFF `
  --config-setting "cmake.define.VCPKG_TARGET_TRIPLET=$arch-windows"
```

Windows wheels place DLLs in `pyalps/bin` alongside the programs; package initialization registers that directory with Python's DLL loader before any extension import. This directory is also included for bindings-only wheels. Downstream Python packages should import `pyalps` before importing their own native extension, so Python can resolve that extension's ALPS DLL dependencies.

## Editable development

After installing the matching SDK, install build dependencies in your active Python environment and use pip's editable mode:

```sh
python -m pip install "scikit-build-core>=1.0" "nanobind>=2.10,<3" "cmake>=4.3" ninja
python -m pip install --no-build-isolation -e python/pyalps
```

Set `ALPS_DIR` as above. On Windows, pass the same `--config-setting` toolchain and triplet arguments as for a wheel build. For a core-only installation, add `--config-setting cmake.define.PYALPS_BUILD_SOLVERS=OFF --config-setting cmake.define.PYALPS_BUNDLE_APPLICATIONS=OFF`.

Python edits take effect in a new interpreter without reinstalling. Rerun the editable install after changing C++ sources, build configuration, or packaged runtime resources. Scikit-build-core owns the Python source mapping; CMake installs the native modules, libraries, XML, downstream headers, and CMake package. Resource lookup uses the installed runtime manifest, so it also works when these files live separately from the Python sources.

## Native runtime layout

The extensions share one nanobind library, named `pyalps_nanobind` to avoid filename collisions with other packages, and one shared ALPS runtime. On Unix, CMake installs relative runtime paths with the package's `lib` directory first and derives any additional dependency directories from resolved link targets. Linux and macOS wheels intended for redistribution must then be repaired with auditwheel or delocate, respectively, to bundle external dependencies and replace build-machine paths; the wheel CI performs this step. See [downstream native extensions](#downstream-native-extensions) for manifest finalization after repair.

## Python compatibility and stable ABI

pyalps ships `cp312-abi3` wheels. Each platform and architecture has one native wheel for GIL-enabled CPython 3.12 and newer; CI builds once and tests that same wheel on Python 3.12, 3.13, and 3.14. Linux libc variants and operating-system deployment targets remain distinct. Python 3.10 and 3.11 are no longer supported.

Every extension and the shared nanobind runtime compile with `STABLE_ABI`, targeting the CPython 3.12 limited API. Packaging CI audits the resulting binaries with abi3audit. Downstream nanobind extensions that exchange ALPS Python objects must also use `STABLE_ABI`; stable-ABI and ordinary nanobind builds have separate type registries. Rebuild existing consumers using the example below.

Free-threaded CPython builds are unsupported: these `abi3` wheels target GIL-enabled CPython, and source builds reject free-threaded interpreters. ALPS relies on the GIL around shared state (`mcobservable`'s reference-count table, the `alps::ngs::signal` singleton, and `mcdata`'s lazily-computed statistics). Supporting free-threading requires making that state thread-safe and selecting the corresponding Python ABI.

## Compatibility and checkpoints

Parameters created from Python retain their Python values. NumPy arrays keep array arithmetic, and changes through a list, array, or shared reference are visible to subsequent Python and C++ reads. A C++ consumer converts the current value to its requested scalar or one-dimensional vector type; incompatible metadata and out-of-range conversions raise an exception. Python metadata may use other shapes and containers supported by the HDF5 writer. Objects such as `None` can be held in memory but have no ALPS HDF5 representation.

The C++ SDK remains independent of Python and nanobind. Python-owned values and their checkpoint decoder are supplied by the bindings. Rebuild downstream C++ extensions against the SDK from the same source revision as the wheel; the parameter layout changed during this migration.

New HDF5 writes distinguish Boolean and signed-byte values with an `__alps_type__` attribute while retaining the existing numeric storage format. Unmarked signed-byte data from old ALPS files retains the legacy Boolean interpretation. The old format cannot distinguish an unmarked `int8` array from a Boolean mask; use a typed reader such as h5py when an old dataset is known to contain signed bytes.

Rectangular mixtures of numeric rows are stored as a single array when every integer remains exact in the common dtype. If mixing integer widths or mixing integers with floating-point or complex rows would round a value, the archive stores the rows separately and reads them back as a list. For example, a `uint64` row containing `2**63 + 1` alongside an `int64` row retains its exact integer values instead of silently converting them to `float64`.

`pyalps.mpi` receives Python objects using matched probes, so asynchronous receives and the wait/test helpers can handle messages larger than mpi4py's default object receive buffer. This adapter exchanges mpi4py messages; Boost.MPI's C++ serialization protocol and skeleton/content API are not wire compatible. Communicating processes must use the same protocol.

## Versioning

pyalps does not carry a version of its own. The numeric version is read from `ALPS_VERSION.txt` at the repository root — the same file `cmake/ALPSVersion.cmake` reads for `ALPS_VERSION_CORE` — so a release bump is one edit rather than two that can drift. `tests/pyalps/test_wheel_payload.py` fails if the installed version and that file disagree.

A prerelease label cannot live in that file: `project(VERSION ...)` rejects a non-numeric version, and neither `find_package()` matching nor the library SOVERSION has a notion of prerelease ordering. CMake takes it from the `ALPS_VERSION_PRERELEASE` cache variable; the Python build takes it from the environment variable of the same name, using the same vocabulary:

| `ALPS_VERSION_PRERELEASE` | version with `ALPS_VERSION.txt` = 2.3.4 |
|---|---|
| unset | `2.3.4` |
| `beta.1` | `2.3.4b1` |
| `alpha.2` | `2.3.4a2` |
| `rc.1` | `2.3.4rc1` |
| `dev.3` | `2.3.4.dev3` |

In GitHub release builds, the provider takes the prerelease label from `GITHUB_REF` (for example, `refs/tags/v3.0.0-beta.1`). It rejects a tag whose numeric version differs from `ALPS_VERSION.txt`, or whose label conflicts with an explicit `ALPS_VERSION_PRERELEASE`. Wheels and source distributions use the same provider. `python python/pyalps/_build_support/alps_version.py` prints the version a build would produce, from any working directory. An sdist preserves its recorded version when rebuilt without the original build environment.

Note the consequence: because the number is inherited, a Python-only API change cannot be signalled in the pyalps version alone — it takes a bump of `ALPS_VERSION.txt`, which moves the whole project.
## Downstream native extensions

The C++ SDK supplies `ALPS::alps` for standalone programs. The installed Python package separately supplies `pyalps::runtime` for extensions that share ALPS objects or HDF5 handles with pyalps. A matching C++ SDK is still required for headers and compile settings.

```cmake
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module Development.SABIModule)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE nanobind_ROOT
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(nanobind 2.10...<3 CONFIG REQUIRED)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pyalps --cmake-dir
  OUTPUT_VARIABLE pyalps_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
find_package(pyalps CONFIG REQUIRED)
nanobind_add_module(my_module NB_STATIC STABLE_ABI my_module.cpp)
target_link_libraries(my_module PRIVATE pyalps::runtime)
```

`pyalps.get_cmake_dir()` exposes the same directory to Python tools. The C++ SDK neither installs this package nor discovers Python. The former SDK function `alps_target_link_pyalps` has been removed.

The same target supplies `<pyalps/export_simulation.hpp>` for exporting a derived simulation through nanobind. This Python-owned header replaces the old SDK header `<alps/ngs/detail/export_sim_to_python.hpp>`; update that include when rebuilding a downstream extension.

Wheel installation writes `pyalps/runtime.json`. After auditwheel or delocate repair, regenerate it with `python python/pyalps/_build_support/runtime_manifest.py --wheel path/to/pyalps.whl`. Cibuildwheel runs this automatically. The manifest records the final relative library paths; pyalps exposes these as imported CMake targets. On macOS, wheel finalization sets linkable `@rpath` library IDs and refreshes their signatures, so downstream builds need no binary-patching commands. Windows uses the matching SDK import libraries and pyalps' registered DLL directory.

CMake derives build-time search paths from the imported targets. If you install or redistribute your extension, set its `INSTALL_RPATH` for the destination layout using normal CMake installation rules. The target does not hard-code the build environment's Python installation into installed extensions. For an extension installed for the same environment, CMake's `INSTALL_RPATH_USE_LINK_PATH` target property can retain the runtime search paths.
