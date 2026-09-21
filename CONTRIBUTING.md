# Contributing to ALPS

Thank you for your interest in ALPS (Algorithms and Libraries for Physics Simulations).
ALPS is a community-driven, open-source ecosystem for numerical simulations of correlated quantum systems.
Contributions at every level — from a one-line bug report to a new simulation method — are welcome and valued.

## Table of contents

- [Ways to contribute](#ways-to-contribute)
- [Reporting bugs and requesting features](#reporting-bugs-and-requesting-features)
- [Getting started with the code](#getting-started-with-the-code)
- [Making a change](#making-a-change)
- [Submitting a pull request](#submitting-a-pull-request)
- [Preparing a release](#preparing-a-release)
- [Review process](#review-process)
- [Code style](#code-style)
- [Recognition](#recognition)
- [Getting help](#getting-help)

---

## Ways to contribute

Contributions fall into four broad levels. You do not need to start at the bottom — jump in wherever your skills fit.

| Level | What this looks like |
|---|---|
| **1 — Feedback** | Install ALPS, try a tutorial, open an issue when something is unclear or broken |
| **2 — Documentation & tutorials** | Improve or extend tutorials on the [ALPS website](https://alps.comp-phys.org), fix documentation errors, add examples |
| **3 — Maintenance** | Fix bugs, improve tests, update dependencies, respond to community questions on Discord |
| **4 — New code** | Contribute a new algorithm, library, or simulation application |

All contributions require agreeing to release your work under the [MIT License](LICENSE.txt).

---

## Reporting bugs and requesting features

Use the [GitHub issue tracker](https://github.com/ALPSim/ALPS/issues). Choose the template that best fits:

- **Bug report** — something is broken or produces wrong results
- **Feature request** — you would like new functionality
- **Simulation help** — you need help setting up a specific model, lattice, or method
- **Website help** — problems with the alps.comp-phys.org website

Before opening a new issue, please search existing issues to avoid duplicates.

---

## Getting started with the code

### Prerequisites

- CMake ≥ 3.27
- A C++17 compiler and C11 compiler (GCC, Clang, or MSVC 2022)
- An installed Boost ≥ 1.76, HDF5 with its C library, and BLAS/LAPACK
- MPI and Boost.MPI for the default parallel build; use `-DALPS_ENABLE_MPI=OFF` for a serial build
- For Fortran examples and simulations: gfortran (or a compatible Fortran compiler).
  The `ALPS_BUILD_FORTRAN` C++ wrapper itself needs no Fortran compiler.
- For Python bindings: Python ≥ 3.10, plus `numpy` and `scipy`

See the [installation page](https://alps.comp-phys.org/install/) for full platform-specific instructions.

### Fork and clone

1. Fork the repository on GitHub.
2. Clone your fork locally:
   ```bash
   git clone https://github.com/<your-username>/ALPS.git
   cd ALPS
   ```
3. Add the upstream remote so you can stay up to date:
   ```bash
   git remote add upstream https://github.com/ALPSim/ALPS.git
   ```

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Alternatively, use the bundled CMake preset:
```bash
cmake --preset default
cmake --build --preset default
```

The Python bindings are a separate `scikit-build-core` project that builds
against an installed ALPS C++ SDK; see the
[`pyalps` build instructions](bindings/python/pyalps/README.md).

Dependencies are discovered through their CMake packages. Set
`CMAKE_PREFIX_PATH` for a non-system installation. ALPS no longer downloads
or compiles a private copy of Boost during configuration, changes the chosen
compiler, or adds `-fpermissive`. The old `Boost_SRC_DIR` and
`ALPS_USE_SYSTEM_BOOST` switches have been removed.

### Native Windows (MSVC)

Install Visual Studio 2022's **Desktop development with C++** workload,
CMake ≥ 3.27, Git, and [vcpkg](https://github.com/microsoft/vcpkg).
Set `VCPKG_ROOT` to its checkout and run in PowerShell:

```powershell
cmake --preset windows-x64
cmake --build --preset windows-x64
ctest --preset windows-x64
cmake --install _build/windows-x64 --config Release
```

The manifest pins Boost, HDF5, OpenBLAS and LAPACK. This preset builds x64
shared libraries and applications with MPI disabled. To build and test Debug,
use the `windows-x64-debug` build and test presets. Install Debug and Release
into separate prefixes (`cmake --install ... --config Debug --prefix ...`).

On a Windows 11 ARM64 host, use `windows-x64-on-arm64` for configuration,
building, and testing (or its `-debug` build/test presets). This opt-in preset
applies a narrowly scoped vcpkg overlay so LAPACK can use x64 GFortran under
Windows emulation. The ordinary x64 preset and CI use upstream ports. The
outputs are x64. Keep Python and all dependencies on the same target
architecture. The install contains the required non-system DLLs in `bin`.

For native Windows ARM64, use `windows-arm64` instead:

```powershell
cmake --preset windows-arm64
cmake --build --preset windows-arm64
ctest --preset windows-arm64
cmake --install _build/windows-arm64 --config Release
```

This preset uses `arm64-windows` dependencies and a
[small numerical-package overlay](cmake/vcpkg-arm64-overlay/README.md) for the
official OpenBLAS ARM64 binaries, including LAPACK 3.12.0. The stock vcpkg
CLAPACK/OpenBLAS combination has incompatible return conventions and fails a
numerical regression. The overlay needs no separate Fortran compiler and uses
the upstream Release C-ABI DLL for both Release and Debug consumers.
Use ARM64 Python for native Python bindings. Keep separate dependency install
directories for x64 and ARM64: vcpkg manifest installation synchronizes its
directory to the requested target and removes packages for other targets.

For a Ninja build, start a matching Visual Studio developer shell and pass
the vcpkg toolchain and triplet explicitly. Build outputs use `bin` for
executables/DLLs and `lib` for link libraries; multi-configuration generators
add their configuration subdirectory automatically.

Keep machine-specific paths, job limits and disk preferences in an untracked
`CMakeUserPresets.json`. To reclaim dependency intermediates automatically, set
`VCPKG_INSTALL_OPTIONS` to
`--clean-buildtrees-after-build;--clean-packages-after-build`. On machines with
limited disk space, setting the Debug executable/shared/module linker flags to
`/DEBUG /INCREMENTAL:NO` retains symbols without large incremental-link caches.

For everyday work, reuse one build directory per configuration and build only
the target being changed, for example:

```powershell
cmake --build --preset windows-x64-debug --target spinmc
```

That builds the target and its dependencies without building every application
and test. The `sdk` preset disables tests, applications and MPI for a small
library build. The default and Windows presets include full native validation.
`BUILD_TESTING` is the single test switch; `ALPS_BUILD_APPLICATIONS` controls
simulation applications and command-line tools together. Examples and tutorial
installation are opt-in. SDK headers are always installed.

`ALPS_BUILD_EXTENSIVE_TESTS=ON` adds the expensive graph and HDF5 type-matrix
tests to `BUILD_TESTING`. The HDF5 matrix compiles each type once and exercises
dataset, attribute and compression modes at runtime; unavailable SZIP encoding
is reported as a skipped test. It replaces `ALPS_BUILD_GRAPH_TESTS` and
`ALPS_BUILD_HDF5_TESTS`.

`ALPS_BUILD_ARCHIVE=ON` adds the optional SQLite archive tool and requires
SQLite3. Its dependency is discovered only when the tool is requested; vcpkg
selects the corresponding manifest feature. The obsolete Boost-source-tree
maintenance tool and its `ALPS_BUILD_DEVELOPER_TOOLS` option have been retired.

`add_subdirectory(ALPS)` defaults to the library alone, with MPI disabled.
An embedding project can explicitly enable the capabilities it needs. MPI,
OpenMP, OpenMP worker scheduling and Fortran remain independent supported
capabilities; worker scheduling requires OpenMP. The unused switch for replacing
the simulation engine's accumulators and the obsolete OpenMPI ULFM prototype
have been retired. The accumulator feature classes used by Python remain.

Migration: replace `ALPS_BUILD_TESTS` with `BUILD_TESTING`, and replace
`ALPS_BUILD_LIBS_ONLY=ON` with `ALPS_BUILD_APPLICATIONS=OFF`. Remove
`ALPS_INSTALL_HEADERS`; select installation components instead if needed.

### Numerical libraries

`BLA_SIZEOF_INTEGER=4` is the default BLAS/LAPACK integer ABI. Use `8` only
with ILP64 dependencies. It replaces the old `LAPACK_64_BIT` alias.
`BLA_VENDOR` and `BLA_STATIC` are passed to CMake's numerical-library
finders. Both numerical libraries are required; missing dependencies cause a
configuration error instead of silently dropping simulation programs. The default build
prefers provider targets to preserve Debug/Release library selection.
The regression suite checks LAPACK's integer ABI and a numerical solve.
`BUILD_TESTING=OFF` also leaves Boost.Test out of the vcpkg manifest features
and the installed SDK never requires it.

### Consuming the C++ SDK

After `cmake --install`, set `CMAKE_PREFIX_PATH` to the SDK prefix (and any
dependency prefixes) or set `ALPS_DIR` to `<prefix>/share/alps`:

```cmake
project(my_simulation LANGUAGES C CXX)
find_package(ALPS CONFIG REQUIRED)
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE ALPS::alps)
```

The exported target carries the include paths, C++17 requirement, compile
definitions and transitive dependencies. Consumers choose their own compiler
and build flags. Use the same ABI and build configuration as the SDK.
Dependency discovery preserves the parent's numerical-provider variables.
`ALPS::headers` exposes the compile interface without linking the library.

The SDK also exports `ALPS::fortran` when the Fortran wrapper is built.
It carries the GNU Fortran compatibility flag needed by the legacy untyped
Fortran bridge; the flag applies only to Fortran consumers of that target.
The two installed Fortran tutorials also require Fortran OpenMP because their
source calls the OpenMP runtime directly.
An SDK with applications exports their executable targets (for example,
`ALPS::spinmc`), listed in `ALPS_APPLICATION_TARGETS`. Consumers may request
`find_package(ALPS CONFIG REQUIRED COMPONENTS applications)` to require them.

Legacy `ALPS_USE_FILE`, `ALPS_LIBRARIES` and dependency-variable aliases have
been removed. Link to the exported targets instead. The C++ package has no
Python discovery or wheel integration. Native Python extensions use the separate
[CMake package supplied by pyalps](bindings/python/pyalps/README.md#downstream-native-extensions).

Installation follows `GNUInstallDirs`, including customized `CMAKE_INSTALL_BINDIR`
and `CMAKE_INSTALL_LIBDIR`. XML resources and optional tutorials live under
`${CMAKE_INSTALL_DATADIR}/alps`, exported as `ALPS_DATA_DIR`.
The former `ALPS_XML_PATH` CMake cache option and `alpsvars` shell scripts have
been removed; the `ALPS_XML_PATH` runtime environment override remains available.

After installing an MPI-disabled LP64 SDK, run the consumer contracts with
`ALPS_DIR=<prefix>/share/alps python -m pytest test/cmake`. They check parent
project defaults, numerical ABI rejection, installed and relocated consumers.
`ALPS_TEST_CMAKE_ARGS` accepts a JSON array of toolchain arguments when needed.

### Run the tests

From the build directory:
```bash
ctest --output-on-failure
```

All tests must pass before submitting a pull request.

---

## Making a change

1. **Sync with upstream** before starting work:
   ```bash
   git fetch upstream
   git checkout master
   git merge upstream/master
   ```

2. **Create a branch** named after what you are doing:
   ```bash
   git checkout -b fix/alea-overflow
   git checkout -b feature/dmrg-excited-states
   git checkout -b docs/tutorial-heisenberg
   ```

3. **Make your changes.** Keep commits focused and self-contained. Write commit messages in the imperative mood:
   ```
   fix: prevent integer overflow in alea accumulator
   feat: add excited-state targeting to DMRG
   docs: add Heisenberg chain tutorial
   ```

4. **Add or update tests** for any changed behaviour. New simulation methods should include at least one regression test comparing output against a known result.

---

## Submitting a pull request

1. Push your branch to your fork:
   ```bash
   git push origin fix/alea-overflow
   ```

2. Open a pull request against the `master` branch of `ALPSim/ALPS`.

3. Fill in the pull request template, including:
   - What problem this solves and why
   - How to test the change
   - Any known limitations or follow-up work

4. Ensure CI passes (build + tests on Linux and macOS).

For substantial changes — new simulation applications, new libraries, significant API modifications — we encourage you to **open an issue or start a discussion first** to get early feedback before investing significant time.

---

## Preparing a release

Update both `ALPS_VERSION.txt` (the C++ SDK version) and `[project].version`
in `pyproject.toml` before creating a release tag. For a final release, both
must be `X.Y.Z` and the tag must be `vX.Y.Z`. For a prerelease such as
`vX.Y.Z-beta.1`, keep the SDK core at `X.Y.Z` and use the Python version
`X.Y.Zb1`. The other supported tag suffixes are `alpha.N`, `rc.N`, and `dev.N`.

Validate the intended tag locally using Python 3.11 or newer:

```bash
python -m pip install packaging
python script/check_release_version.py --ref refs/tags/vX.Y.Z
```

The packaging workflow checks these versions before building and checks every
wheel and source distribution, including its embedded metadata, before upload.
Tag pushes publish the full release to PyPI, including CPython 3.9–3.14 wheels.
Merge and validate the release commit before tagging it. Keep tags fixed once
their release has been published.

If a published tag contains the wrong version, rerunning its workflow will
rebuild the same incorrect artifacts. Correct both version files first. If
the intended version has no distributions on PyPI, maintainers can approve
resetting the tag to the validated correction and publishing that version.
If the intended version already has distributions, prepare a new patch
release instead: PyPI does not allow replacing uploaded filenames. Do not
use `skip-existing` to hide a version mismatch.

---

## Review process

ALPS uses a consensus-based review model:

- Pull requests are reviewed by **maintainers** (at least one per simulation code) and **core maintainers**.
- A pull request is accepted if all active reviewers approve, or if no objections are raised within **six weeks** of submission.
- Controversial changes can be escalated to the [Governing Council](https://alps.comp-phys.org/govern/).

Core maintainers are responsible for validating that code compiles, tests pass, and results are physically correct. Please be responsive to review comments; PRs with no author activity for eight weeks may be closed.

If you are contributing a new simulation application or library, the Governing Council will discuss a maintenance commitment with you — typically a few hours per month for bug fixes, dependency updates, and community support.

---

## Code style

### C++

- Target C++17.
- Match the style of the surrounding code. ALPS does not enforce a single formatter, but keeps consistent conventions within each subdirectory.
- Avoid undefined behaviour and compiler warnings. New code should compile cleanly with `-Wall -Wextra` on GCC and Clang.
- Prefer standard library and Boost facilities over hand-rolled implementations.

### Python

- Follow [PEP 8](https://peps.python.org/pep-0008/).
- Type annotations are encouraged for new public functions.

### CMake

- CMake ≥ 3.27 features are acceptable. Express dependencies and compiler settings
  on targets with explicit `PRIVATE`, `PUBLIC` or `INTERFACE` scope.
- Use target-based linking (`target_link_libraries`, `target_include_directories`) rather than directory-level commands.

---

## Recognition

ALPS releases are accompanied by a publication in a peer-reviewed journal. **Active contributors are added as co-authors.** The Governing Council decides the author list for each release, taking into account contributions to code, documentation, tutorials, testing, and community support.

Contributing documentation, tutorials, or code (Level 2 — improving or extending tutorials and website documentation — or above) with sustained effort is the typical threshold for co-authorship consideration.

---

## Getting help

| Channel | Use it for |
|---|---|
| [Discord](https://discord.gg/JRNWnnva9g) | Questions about using ALPS, development discussion, meeting the community |
| [GitHub Issues](https://github.com/ALPSim/ALPS/issues) | Bug reports, feature requests, concrete problems with the code |
| [ALPS website](https://alps.comp-phys.org) | Documentation, tutorials, governance, events |
| [Governing Council](https://alps.comp-phys.org/govern/) | Onboarding for new simulation codes, co-authorship, major contributions |

We look forward to your contribution!
