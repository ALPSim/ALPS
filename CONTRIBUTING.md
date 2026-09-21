# Contributing to ALPS

Thank you for your interest in ALPS (Algorithms and Libraries for Physics Simulations). ALPS is a community-driven, open-source ecosystem for numerical simulations of correlated quantum systems. Contributions at every level — from a one-line bug report to a new simulation method — are welcome and valued.

For release history and migration guidance, see [CHANGELOG.md](CHANGELOG.md).

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

- CMake ≥ 4.3; Ninja for the `default`, `sdk`, and `distribution` presets ([installation instructions](#install-cmake-and-ninja))
- A C++17 compiler and C11 compiler (GCC, Clang, or MSVC 2022)
- An installed Boost ≥ 1.76, HDF5 with its C library, and BLAS/LAPACK
- MPI and Boost.MPI when configuring with `-DALPS_ENABLE_MPI=ON`
- For Fortran examples and simulations: gfortran (or a compatible Fortran compiler). The SDK always includes the C++ bridge `ALPS::fortran`; building it needs no Fortran compiler.
- For Python bindings: GIL-enabled CPython ≥ 3.12, plus `numpy` and `scipy`

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

### Install CMake and Ninja

Check `cmake --version` first: ALPS requires **4.3 or newer**. Your system package manager or IDE may provide an older version. If you need an upgrade, installing [CMake through pip](https://cmake-python-distributions.readthedocs.io/en/stable/installation.html) in a virtual environment is a convenient option on Linux, macOS, and Windows. The [CMake wheels](https://pypi.org/project/cmake/#files) include binaries for x64 and ARM64 on all three platforms.

With Python 3.12 or newer installed, run the following from the repository root. If you already have an active virtual environment, skip creating a new one and run the two pip commands in it.

Linux and macOS (bash or zsh):

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --upgrade "cmake>=4.3" ninja
```

If Linux reports that `venv` or `ensurepip` is unavailable, install your distribution's Python venv package first (for example, `python3-venv` on Debian/Ubuntu).

Windows (PowerShell):

```powershell
py -3 -m venv .venv
$env:PATH = "$(Resolve-Path .venv\Scripts);$env:PATH"
python -m pip install --upgrade pip
python -m pip install --upgrade "cmake>=4.3" ninja
```

The PowerShell PATH assignment makes the environment's tools available in the current shell without running an activation script. In each new terminal, repeat the activation command on Linux/macOS or the PATH assignment on Windows before building. Configure your IDE to use the same CMake executable if needed.

Verify the tools in the shell you will use for the build:

```sh
cmake --version
ctest --version
ninja --version
```

Both CMake and CTest must report 4.3 or newer. If an older executable still wins, check its location with `command -v cmake` on Linux/macOS or `Get-Command cmake` in PowerShell and put the virtual environment's executable directory first on `PATH`.

You can also use [official CMake downloads](https://cmake.org/download/) without Python: install or unpack a release ≥ 4.3 for your OS and architecture, and add its executable directory to `PATH`. For the macOS application bundle this is `/Applications/CMake.app/Contents/bin`. Install Ninja separately if using a Ninja preset; the Windows Visual Studio presets do not require it. A system package is equally suitable when `cmake --version` confirms it meets the requirement.

### Build

Shared C++ libraries live under `src/alps/`; command-line utilities live under `src/tools/`. Solver implementations, including MaxEnt, currently live under `src/apps/`. Python bindings have their own project under `python/pyalps/`. Bundled dependencies live under `third_party/`: the temporary [Numeric Bindings headers](third_party/boost_numeric_bindings/README.md) and the [XDR serialization implementation](third_party/xdr/README.md). Shared XML definitions and stylesheets live under `src/alps/resources/`, tests under `tests/`, and CI/release helpers under `.github/scripts/`. The [tutorial guide](tutorials/README.md) is the learning entry point; solver references live beside the tutorials, and focused library examples live under `tutorials/examples/`.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

CTest supplies the source XML resources automatically. When running an uninstalled example or executable directly, set `ALPS_XML_PATH` to the absolute path of `src/alps/resources/`. Installed programs use the SDK's `share/alps/xml/` resources.

Alternatively, use the bundled CMake preset:
```bash
cmake --preset default
cmake --build --preset default
```

The Python bindings are a separate `scikit-build-core` project that builds against an installed ALPS C++ SDK; see the [`pyalps` build instructions](python/pyalps/README.md).

Tutorial sources are an opt-in installation component. After installing the SDK with `cmake --install build`, use `cmake --install build --component tutorials` to add them under `share/alps/tutorials`; no configuration switch is needed. This includes the library examples and their standalone SDK build instructions. The Python extension example belongs to the [`pyalps` package](python/pyalps/examples/ising/README.md).

Dependencies are discovered through their CMake packages. Set `CMAKE_PREFIX_PATH` for a non-system installation. The presets select Release; plain CMake invocations use the configuration you specify. Language requirements propagate through the SDK targets, so consumers can select C++17 or a newer standard.

### Native Windows (MSVC)

Install Visual Studio 2022's **Desktop development with C++** workload, [CMake ≥ 4.3](#install-cmake-and-ninja), Git, and [vcpkg](https://github.com/microsoft/vcpkg). Set `VCPKG_ROOT` to its checkout. Use the preset matching your native architecture. On Windows x64, run in PowerShell:

```powershell
cmake --preset windows-x64
cmake --build --preset windows-x64
ctest --preset windows-x64
cmake --install _build/windows-x64 --config Release
```

The manifest pins Boost, HDF5, OpenBLAS and LAPACK. This preset builds x64 shared libraries and applications with MPI disabled. To build and test Debug, pass `--config Debug` to the build command and `-C Debug` to CTest. Install Debug and Release into separate prefixes (`cmake --install ... --config Debug --prefix ...`).

On Windows ARM64, use `windows-arm64`:

```powershell
cmake --preset windows-arm64
cmake --build --preset windows-arm64
ctest --preset windows-arm64
cmake --install _build/windows-arm64 --config Release
```

This preset uses `arm64-windows` dependencies and a [small numerical-package overlay](cmake/vcpkg-arm64-overlay/README.md) for the official OpenBLAS ARM64 binaries, including LAPACK 3.12.0. The overlay supplies compatible BLAS/LAPACK interfaces, needs no separate Fortran compiler, and uses the upstream Release C-ABI DLL for both Release and Debug consumers. Use Python and dependencies matching the target architecture: x64 for `windows-x64`, ARM64 for `windows-arm64`. The install contains the required non-system DLLs in `bin`. Keep separate dependency install directories for x64 and ARM64: vcpkg manifest installation synchronizes its directory to the requested target and removes packages for other targets.

For a Ninja build, start a matching Visual Studio developer shell and pass the vcpkg toolchain and triplet explicitly. Build outputs use `bin` for executables/DLLs and `lib` for link libraries; multi-configuration generators add their configuration subdirectory automatically.

Keep machine-specific paths, job limits and disk preferences in an untracked `CMakeUserPresets.json`. To reclaim dependency intermediates automatically, set `VCPKG_INSTALL_OPTIONS` to `--clean-buildtrees-after-build;--clean-packages-after-build`. On machines with limited disk space, setting the Debug executable/shared/module linker flags to `/DEBUG /INCREMENTAL:NO` retains symbols without large incremental-link caches.

For everyday work, reuse one build directory per configuration and build only the target being changed, for example:

```powershell
cmake --build --preset windows-x64 --config Debug --target spinmc
```

That builds the target and its dependencies without building every application and test. The `sdk` preset disables tests, applications and MPI for a small library build. The default and Windows presets include full native validation. `BUILD_TESTING` is the single test switch; `ALPS_BUILD_APPLICATIONS` controls simulation applications and command-line tools together. Examples and tutorial installation are opt-in. SDK headers are always installed.

`ALPS_BUILD_EXTENSIVE_TESTS=ON` adds the expensive graph and HDF5 type-matrix tests to `BUILD_TESTING`. The HDF5 matrix compiles each type once and exercises dataset, attribute and compression modes at runtime; unavailable SZIP encoding is reported as a skipped test.

`add_subdirectory(ALPS)` defaults to the library alone, with MPI disabled. An embedding project can explicitly enable the capabilities it needs. MPI is opt-in in every build. `ALPS_ENABLE_OPENMP=ON` enables OpenMP, including worker scheduling; select the worker's thread count at runtime.

`ALPS_BUILD_EXAMPLES=ON` builds the C++ examples. The Fortran examples are a separate consumer of the installed SDK, so ordinary builds need no Fortran compiler:

```bash
cmake -S tutorials/examples/fortran -B build/fortran -DCMAKE_PREFIX_PATH=/path/to/alps
cmake --build build/fortran
ctest --test-dir build/fortran --output-on-failure
```

### Numerical libraries

The SDK uses LP64 BLAS/LAPACK: 32-bit integers and lowercase symbols with a trailing underscore. Alternate integer widths and symbol spellings are unsupported. `BLA_VENDOR` and `BLA_STATIC` are passed to CMake's numerical-library finders. Both numerical libraries are required; missing dependencies cause a configuration error. The default build prefers provider targets to preserve Debug/Release library selection. The regression suite checks LAPACK's integer ABI and a numerical solve. `BUILD_TESTING=OFF` also leaves Boost.Test out of the vcpkg manifest features and the installed SDK never requires it.

### Consuming the C++ SDK

After `cmake --install`, set `CMAKE_PREFIX_PATH` to the SDK prefix (and any dependency prefixes) or set `ALPS_DIR` to `<prefix>/share/alps`:

```cmake
project(my_simulation LANGUAGES C CXX)
find_package(ALPS CONFIG REQUIRED)
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE ALPS::alps)
```

The exported target carries the include paths, C++17 requirement, compile definitions and transitive dependencies. Consumers choose their own compiler and build flags. Use the same ABI and build configuration as the SDK. Dependency discovery preserves the parent's numerical-provider variables. `ALPS::headers` exposes the compile interface without linking the library.

The SDK also exports the C++ Fortran bridge as `ALPS::fortran`. It carries the GNU Fortran compatibility flag needed by the legacy untyped Fortran bridge; the flag applies only to Fortran consumers of that target. The two installed Fortran tutorials also require Fortran OpenMP because their source calls the OpenMP runtime directly. An SDK with applications exports their executable targets (for example, `ALPS::spinmc`), listed in `ALPS_APPLICATION_TARGETS`. Consumers may request `find_package(ALPS CONFIG REQUIRED COMPONENTS applications)` to require them.

The same build installs the static solver libraries `ALPS::maxent`, `ALPS::cthyb`, and `ALPS::ctint`. Require `COMPONENTS solvers` and include `<alps/solvers.hpp>` to run them from C++. Each function accepts `alps::params` and an output filename. The native programs and Python wrappers use these same implementations; building pyalps never compiles application sources.

Native Python extensions use the separate [CMake package supplied by pyalps](python/pyalps/README.md#downstream-native-extensions).

Installation follows `GNUInstallDirs`, including customized `CMAKE_INSTALL_BINDIR` and `CMAKE_INSTALL_LIBDIR`. XML resources and optional tutorials live under `${CMAKE_INSTALL_DATADIR}/alps`, exported as `ALPS_DATA_DIR`. Use the `ALPS_XML_PATH` runtime environment override to select different XML resources.

After installing an MPI-disabled LP64 SDK, run the consumer contracts with `ALPS_DIR=<prefix>/share/alps python -m pytest tests/cmake`. They check parent project defaults, numerical ABI rejection, installed and relocated consumers. `ALPS_TEST_CMAKE_ARGS` accepts a JSON array of toolchain arguments when needed.

### XML command-line tools

Unix installations provide `alps-xml` for rendering simulation results and extracting plots. It requires Python 3 and `xsltproc` on `PATH` (on Ubuntu, `sudo apt install xsltproc`; on macOS, `brew install libxslt` and add its `bin` directory to `PATH`). It finds the installed stylesheets relative to the executable, including after the installation is moved.

```bash
alps-xml plot text results.plot.xml
alps-xml plot grace results.plot.xml --output results.agr
alps-xml convert html simulation.out.xml --output results.html
alps-xml extract text plot-definition.xml task*.out.xml --output measurements.txt
```

Plot and extraction formats are `text`, `html`, `gnuplot`, `matplotlib` and `grace`. Conversion formats are `text` and `html`. Extraction accepts simulation files or XML archives. The generated Matplotlib program runs with Python 3. The `xml` install component contains both this command and its resources.

### Run the tests

From the build directory:
```bash
ctest --output-on-failure
```

All tests must pass before submitting a pull request.

---

The XML CLI integration tests install the `xml` component into a temporary prefix, relocate it and exercise real transformations. Run them against a configured Unix build with applications enabled and `xsltproc` available:

```bash
ALPS_XML_BUILD=build python -m pytest tests/cli -q
```

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

5. **Update the Unreleased section of [CHANGELOG.md](CHANGELOG.md)** for user-facing features, fixes, removals, build requirements, or migration steps. Internal changes without a user-facing effect do not need an entry.

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

4. Ensure CI passes, including source builds, native Windows builds, and Python wheel validation.

For substantial changes — new simulation applications, new libraries, significant API modifications — we encourage you to **open an issue or start a discussion first** to get early feedback before investing significant time.

---

## CI coverage

Pull requests against any branch, merge-queue entries, and pushes to `master` run four source configurations: GCC 11 with the minimum Boost 1.76, GCC 14 with installed examples and XML tools, Clang 18 with C++23, and Apple Clang on macOS ARM64. The separate Windows workflow retains x64 Release/Debug and ARM64 Release builds. Packaging workflows test the same `cp312-abi3` wheel across Python 3.12–3.14, including native Windows and macOS forward compatibility.

The source workflow runs its full fourteen-configuration matrix weekly and on release tags. The full tier adds GCC 15, Clang 14 and 22, Linux ARM64, macOS Intel and macOS 26, representative intermediate Boost releases, C++20 with the extensive graph and HDF5 tests, and AddressSanitizer plus UndefinedBehaviorSanitizer. Sanitizer runs disable MPI and dependency leak detection; address and undefined-behavior errors fail the job. This tests supported boundaries and representative combinations without rebuilding every compiler/Boost permutation on every pull request.

To request the full matrix before merging, open **Actions → ALPS source CI → Run workflow**, select the branch and the `full` tier. The single manifest [`.github/ci-matrix.json`](.github/ci-matrix.json) defines both tiers and pins the Boost archive checksums. The `Source CI` check aggregates matrix results and is suitable as a required branch check. Workflow linting and CI-helper tests run before compilation; native and installed-wheel jobs retain test reports and display result counts in their job summaries. Manual packaging runs build and test artifacts; only a pushed release tag can publish to PyPI.

## Preparing a release

Review the Unreleased entries in [CHANGELOG.md](CHANGELOG.md), group related changes, and check the migration guidance. At release time, give the section the release version and date, then start a new Unreleased section for subsequent work.

Update `ALPS_VERSION.txt`, the shared SDK and Python release version, before creating a release tag. For a final release, it must be `X.Y.Z` and the tag must be `vX.Y.Z`. For a prerelease such as `vX.Y.Z-beta.1`, keep the file at `X.Y.Z`; the Python metadata provider derives `X.Y.Zb1` from the tag. The other supported tag suffixes are `alpha.N`, `rc.N`, and `dev.N`.

Validate the intended tag locally using Python 3.12 or newer:

```bash
python -m pip install packaging
python .github/scripts/check_release_version.py --ref refs/tags/vX.Y.Z
```

The packaging workflow checks these versions before building and checks every wheel and source distribution, including its embedded metadata, before upload. Tag pushes publish the full release to PyPI, including `cp312-abi3` wheels shared by GIL-enabled CPython 3.12 and newer. CI tests the same wheel on 3.12, 3.13, and 3.14. Merge and validate the release commit before tagging it. Keep tags fixed once their release has been published.

If a published tag contains the wrong version, rerunning its workflow will rebuild the same incorrect artifacts. Correct `ALPS_VERSION.txt` and validate the tag first. If the intended version has no distributions on PyPI, maintainers can approve resetting the tag to the validated correction and publishing that version. If the intended version already has distributions, prepare a new patch release instead: PyPI does not allow replacing uploaded filenames. Do not use `skip-existing` to hide a version mismatch.

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

- CMake ≥ 4.3 features are acceptable. Express dependencies and compiler settings on targets with explicit `PRIVATE`, `PUBLIC` or `INTERFACE` scope.
- Use target-based linking (`target_link_libraries`, `target_include_directories`) rather than directory-level commands.

### Markdown

Keep each prose paragraph or list item's text on one source line and let the renderer wrap it. Preserve blank lines between paragraphs and the line structure of code blocks, tables, and nested lists.

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
