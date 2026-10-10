# Contributing to ALPS

Bug reports, documentation, tests and new simulation methods are welcome. To use ALPS without developing it, start with the [installation instructions](README.md#installation) and [tutorial guide](tutorials/README.md). See [CHANGELOG.md](CHANGELOG.md) for unreleased changes and migration guidance.

- [Getting started with the code](#getting-started-with-the-code)
- [Run the tests](#run-the-tests)
- [Making a change](#making-a-change)
- [Provenance and scientific credit](#provenance-and-scientific-credit)
- [Build reference](#build-reference)
- [CI coverage](#ci-coverage)
- [Preparing a release](#preparing-a-release)

## Reporting bugs and requesting features

Search the [issue tracker](https://github.com/ALPSim/ALPS/issues) before opening an issue. For a bug, include your ALPS commit or version, OS and architecture, compiler and dependency versions, exact commands, the error output, and a small reproducer. For incorrect simulation results, include the input parameters, expected result and its reference, and any relevant random seed.

Small fixes can go directly to a pull request. Discuss a substantial new application, library or API change with the maintainers before implementing it; see the project's [contribution policy](https://alps.comp-phys.org/govern/contribute/). Contributions must be compatible with the [MIT License](LICENSE.txt).

## Getting started with the code

Build and install the C++ SDK with CMake, then build the Python package against it with pip. The workflow below includes applications and tests, uses one SDK installation, and requires no environment manager. Run commands from the repository root unless stated otherwise.

### Fork and clone

Fork [ALPSim/ALPS](https://github.com/ALPSim/ALPS) on GitHub, replace `<your-username>` below, and clone your fork:

```sh
git clone https://github.com/<your-username>/ALPS.git
cd ALPS
git remote add upstream https://github.com/ALPSim/ALPS.git
```

### Prerequisites

- CMake ≥ 3.27, Ninja for the bundled presets, and C++17/C11 compilers such as GCC or Clang.
- Boost ≥ 1.76 with its compiled libraries and CMake packages, HDF5 ≥ 1.10.5 (C library), and LP64 BLAS/LAPACK. Use serial HDF5 for the default MPI-disabled build; see [numerical libraries](#numerical-libraries).
- For Python development: GIL-enabled CPython ≥ 3.11 in a writable Python environment. Pip installs NumPy, SciPy and Matplotlib with pyalps. Free-threaded Python is unsupported.
- Optional: MPI and Boost.MPI for `ALPS_ENABLE_MPI=ON`; an OpenMP runtime for `ALPS_ENABLE_OPENMP=ON`; a Fortran compiler for the Fortran examples.

Use existing dependencies or your preferred package manager. Keep the compiler, architecture and native dependency stack consistent between the SDK, bindings and downstream extensions. Older website instructions for a combined Boost.Python build do not describe this checkout.

For example, on Ubuntu 24.04:

```sh
sudo apt-get update
sudo apt-get install build-essential libboost-all-dev libhdf5-dev libblas-dev liblapack-dev
```

On macOS, install Apple's Command Line Tools with `xcode-select --install` if needed. If you use Homebrew:

```sh
brew install boost hdf5 openblas
export CMAKE_PREFIX_PATH="$(brew --prefix boost):$(brew --prefix hdf5):$(brew --prefix openblas)"
```

These are optional provider examples. For another non-system installation, set the `CMAKE_PREFIX_PATH` environment variable to its dependency prefixes, separated by colons on Linux/macOS. Keep it set for both SDK and Python builds. The CMake command-line form instead uses semicolons: `-DCMAKE_PREFIX_PATH="/prefix/one;/prefix/two"`.

### Windows users

On Windows, use a Linux distribution inside WSL and follow the Linux instructions here.

### Install CMake and Ninja

Check `cmake --version`, `ctest --version` and `ninja --version`. Reuse suitable tools and an existing Python environment. If you need a Python environment and build tools, one option is:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "cmake>=3.27" ninja numpy h5py
```

Skip the first two lines if already using a suitable environment. On Debian/Ubuntu, creating a venv may first require `sudo apt-get install python3-venv`. Activate the same environment in each new terminal. Both CMake and CTest must be at least 3.27; use `command -v cmake` to check which installation your shell finds.

Application tests require a Python ≥ 3.11 interpreter; MaxEnt reference tests also require NumPy and h5py. CMake and Ninja can also come from system packages or [official CMake binaries](https://cmake.org/download/). A generator other than Ninja can be selected with a plain CMake invocation instead of a preset.

### Build

Citation metadata and application rules live in `CITATION.cff` and `CITATIONS.yaml`.
See [citation maintenance](.github/scripts/citations/README.md) for generation and validation.
Native builds use checked-in generated citation data and do not require Python.
Editing that data requires Python ≥ 3.9 with `PyYAML` and `jsonschema`:
`python -m pip install -r .github/scripts/citations/requirements.txt`, then
`python .github/scripts/generate_citations.py --regenerate`. Commit the generated files
alongside the authorities; CI checks that they agree. For the optional citation
tests, select an interpreter with `-DALPS_CITATION_PYTHON=/path/to/python`.


Build, test and install the shared SDK and applications:

```sh
cmake --preset default
cmake --build --preset default --parallel 2
ctest --preset default
cmake --install _build/default
export ALPS_DIR="$PWD/_build/default/install/share/alps"
export PATH="$PWD/_build/default/install/bin:$PATH"
```

The preset selects Release and installs into `_build/default/install`, without administrator privileges. `ALPS_DIR` selects that SDK for Python and downstream CMake builds; it does not replace the dependency prefixes above. Keep these exports in each development shell. Two parallel compile jobs are a conservative starting point; adjust to your available memory.

For Python development, continue in your active Python environment:

```sh
python -m pip install "scikit-build-core>=1.0" "nanobind==2.15.0" \
  "patchelf>=0.14; sys_platform == 'linux'"
CMAKE_BUILD_PARALLEL_LEVEL=2 python -m pip install --no-build-isolation -e python/pyalps \
  --config-setting "build-dir=$PWD/_build/python"
python -c "import pyalps; print(pyalps.__file__)"
```

This editable installation reuses `_build/python` for binding builds. `--no-build-isolation` requires the build dependencies installed above; the nanobind pin must match the package's [build requirements](python/pyalps/pyproject.toml). The default bindings include solver modules and bundled programs, so they require a shared SDK with applications enabled, as built here. The [Python package guide](python/pyalps/README.md) covers wheels and smaller core-only builds; its `distribution` preset is an alternative SDK build, not an additional prerequisite for this workflow.

### Edit, rebuild and rerun

- **Python source:** edits take effect in a new interpreter without reinstalling.
- **SDK or application C++:** rerun the SDK build and install commands, then the editable pip command if using Python. The Python installation contains copies of native runtime files.
- **Binding C++, package resources or build configuration:** rerun the editable pip command.
- **Compiler, SDK or dependency-provider change:** reconfigure the SDK deliberately and add `--config-setting cmake.args=--fresh` to the next editable install to clear stale binding discovery. Do not reuse objects from an incompatible toolchain.

Reuse build directories for ordinary edits. Machine-specific CMake settings belong in an untracked `CMakeUserPresets.json`. For a focused native rebuild, use `cmake --build --preset default --target <target> --parallel 2`.

## Run the tests

Native MaxEnt reference tests require NumPy and h5py in the CMake-selected Python interpreter (`python -m pip install numpy h5py`). The executable and Python binding share the same scientific validation.

CTest runs the native suite only. After the SDK build, run:

```sh
ctest --preset default
```

After the editable install above, test Python and installed-SDK consumers too:

```sh
python -m pip install "pytest>=8"
PYALPS_TEST_DOWNSTREAM_EXPORT=1 CMAKE_BUILD_PARALLEL_LEVEL=2 \
  python -m pytest python/pyalps/tests tests/cmake -q -rs
```

Keep `ALPS_DIR` and any dependency prefixes set. The downstream flag enables tests that compile Python extensions against the installed SDK and pyalps runtime. These checks compile additional small projects and take longer than import tests. `tests/cmake` assumes an MPI-disabled LP64 SDK for its consumer contracts; use the default SDK for this command. If custom toolchain arguments are needed by these temporary builds, `ALPS_TEST_CMAKE_ARGS` accepts a JSON array of CMake arguments.

Read the skip reasons: MPI tests need additional MPI setup, and some wheel checks apply only to repaired distribution artifacts. A successful local run with skips does not exercise every CI configuration. For a quick iteration, select native tests with `ctest --preset default -R <pattern>` or Python tests with `python -m pytest <test-file> -q`.

For XML CLI changes, install `xsltproc` (Ubuntu: `sudo apt-get install xsltproc`; Homebrew: `brew install libxslt` and add its `bin` directory to `PATH`), then run:

```sh
ALPS_XML_BUILD="$PWD/_build/default" python -m pytest tests/cli -q
```

These tests install and relocate the XML component before exercising transformations. Build/release helper changes also have tests under `tests/ci` and `tests/packaging`; release-version helpers require Python ≥ 3.11 and `packaging`. Run the tests relevant to your change, report failures or skipped coverage, and let CI validate the broader platform matrix.

## Making a change

Start from an up-to-date local `master`, then create one topic branch:

```sh
git fetch upstream
git switch master
git merge --ff-only upstream/master
git switch -c fix/alea-overflow
```

Replace `fix/alea-overflow` with your own branch name. If the fast-forward fails, resolve your local branch history before continuing; do not discard unrelated work.

Keep commits focused. Add or update tests for changed behavior; scientific changes should include a small regression against a known result. Update documentation and the Unreleased section of [CHANGELOG.md](CHANGELOG.md) for user-facing features, fixes, removals or migration steps. Internal changes without a user-facing effect do not need changelog entries.

### Code style

Match the surrounding code. Target C++17, avoid undefined behavior, and check compiler warnings. Follow [PEP 8](https://peps.python.org/pep-0008/) for Python. CMake changes must work with 3.27 and express usage requirements on targets with explicit `PRIVATE`, `PUBLIC` or `INTERFACE` scope. Let the CMake build configuration or caller select optimization flags; examples and tests must not impose custom optimization settings. In Markdown, keep each prose paragraph on one source line while preserving code blocks, tables and list structure.

### Submitting a pull request

Commit your changes and push the current topic branch to your fork:

```sh
git push -u origin HEAD
```

Open a pull request against `ALPSim/ALPS:master`. Explain the problem, resulting behavior, tests run and any limitations using the PR template. Respond to review comments and address relevant CI failures. Review, maintenance commitments and contributor recognition follow the published [contribution policy](https://alps.comp-phys.org/govern/contribute/) and [governance](https://alps.comp-phys.org/govern/).

For development questions, use [Discord](https://discord.gg/JRNWnnva9g); reproducible bugs belong in the [issue tracker](https://github.com/ALPSim/ALPS/issues).

## Provenance and scientific credit

These expectations apply to human and AI-assisted contributions alike. Record
provenance while making the change, when the sources are known.

- When copying, translating, or substantially adapting external code, add a
  comment near the affected code identifying the upstream project, source file,
  and version or commit where available. Describe the relationship accurately
  (for example, copied, translated, or adapted).
- Preserve existing copyright and license notices, and include any required
  upstream license text with third-party material. Identify that material and
  its terms in the pull request for maintainer review. Flag uncertain provenance
  or licensing before merge; do not assume that ALPS's MIT license replaces
  upstream terms.
- Credit the original method papers and upstream implementations that a new or
  changed component builds on. Update bibliographic records in
  [CITATION.cff](CITATION.cff) and the relevant component mappings in
  [CITATIONS.yaml](CITATIONS.yaml) in the same pull request, following the
  [citation maintenance instructions](.github/scripts/citations/README.md). Scientific
  credit is separate from license compliance; references should be relevant to
  the affected component.
- Do not invent attribution or claim independent implementation without
  evidence. State what is known and flag gaps for review.

Maintainers review provenance and citation changes as part of normal pull
request review.

---

## Build reference

### Repository layout

- `src/alps/`: C++ components, public headers and module-local tests; see the [module map and library boundaries](src/alps/README.md).
- `src/apps/` and `src/tools/`: simulation applications, shared solver implementations and [command-line tools](src/tools/README.md).
- `python/pyalps/`: Python sources, bindings, packaging and extension example.
- `tutorials/`: [ordered tutorials and standalone library examples](tutorials/README.md).
- `tests/`: cross-module integration, Python, SDK, CLI, build-helper and packaging tests.
- `third_party/`: [Numeric Bindings headers](third_party/boost_numeric_bindings/README.md) and [XDR serialization](third_party/xdr/README.md).
- `cmake/` and `.github/`: shared build configuration, generated-header templates, version file, CI and release helpers.

Add exported headers to the owning target's CMake `HEADERS` file set. Public `<alps/...>` and `<ietl/...>` include names are independent of the physical source directory. Keep subsystem tests beside their owner and cross-module tests under `tests/`.

### Build options

| Option | Default for a top-level build | Purpose |
| --- | --- | --- |
| `ALPS_BUILD_TESTING` | `ON` | Build and register native ALPS tests, independently of a parent's `BUILD_TESTING` |
| `ALPS_BUILD_APPLICATIONS` | `ON` | Build simulation applications, solver libraries and command-line tools |
| `BUILD_SHARED_LIBS` | `ON` when unset | Shared libraries; required for Python bindings |
| `ALPS_ENABLE_MPI` | `OFF` | Enable MPI; requires matching MPI and Boost.MPI installations |
| `ALPS_ENABLE_OPENMP` | `OFF` | Enable OpenMP, including worker scheduling |
| `ALPS_BUILD_EXTENSIVE_TESTS` | `OFF` | Add expensive graph and HDF5 type-matrix tests when testing is enabled |

For example, configure with `cmake --preset default -DALPS_ENABLE_OPENMP=ON`. The `sdk` preset disables applications and tests; `distribution` disables tests and its build preset installs automatically. When embedding ALPS with `add_subdirectory`, applications and tests default to `OFF`. MPI remains opt-in. Headers and the C++ Fortran bridge are always part of the SDK; building that bridge needs no Fortran compiler.

Examples build separately against an installed SDK; see the [C++ and Fortran example instructions](tutorials/examples/README.md). Fortran tutorials that call OpenMP also need a Fortran OpenMP runtime. To install tutorial sources under `share/alps/tutorials`, run `cmake --install _build/default --component tutorials` after installing the SDK.

### Numerical libraries

Both BLAS and LAPACK are required, with LP64 (32-bit) integers and lowercase symbols ending in an underscore. ILP64 and alternate symbol spellings are unsupported. `BLA_VENDOR` and `BLA_STATIC` select a provider or static numerical libraries through CMake's finders. Keep this ABI consistent when building downstream consumers; installing the SDK does not supply the external numerical libraries.

### Consuming the C++ SDK

Use the installed `ALPS_DIR` from the build instructions, or add the SDK installation prefix to `CMAKE_PREFIX_PATH`, alongside dependency prefixes:

```cmake
cmake_minimum_required(VERSION 3.27)
project(my_simulation LANGUAGES C CXX)
find_package(ALPS CONFIG REQUIRED)
add_executable(my_simulation main.cpp)
target_link_libraries(my_simulation PRIVATE ALPS::alps)
```

The imported target carries headers, C++17 requirements, compile definitions and transitive dependencies. Use a compiler and configuration compatible with the SDK's ABI. `ALPS::headers` exposes the compile interface without linking; `ALPS::fortran` supplies the C++ Fortran bridge and its GNU Fortran compatibility flag. Programs using only utilities can request `find_package(ALPS CONFIG REQUIRED COMPONENTS utilities)` and link `ALPS::utilities`. Archive-only programs can similarly request the `hdf5` component and link `ALPS::hdf5`, which also owns archive signal cleanup. Typed parameter programs can request the `params` component and link `ALPS::params`. XML parsing/output and command-line parsing are available through the `xml` and `cli` components and targets `ALPS::xml` and `ALPS::cli`. The parameter-file constructor and text/XML conversion adapters still require `ALPS::alps`. Package discovery still checks the SDK's complete dependency set; component-specific configuration is future work.

Use `ALPS::containers` for array storage and `ALPS::numerics` for matrix/vector algorithms and array mathematics, including the existing `<alps/multi_array.hpp>` umbrella. Numerical archive consumers link `ALPS::numeric_io` and explicitly include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>`. For example:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS numeric_io)
target_link_libraries(my_simulation PRIVATE ALPS::numeric_io)
```

`<alps/numeric/matrix.hpp>` no longer includes its HDF5 adapter automatically. See the [migration notes](CHANGELOG.md#migration) for this include requirement and the `<alps/xml/matrix.hpp>` adapter retaining matrix XML output.

An SDK built with applications also exports executable targets such as `ALPS::spinmc` and the solver libraries `ALPS::maxent`, `ALPS::cthyb` and `ALPS::ctint`. Require them with `find_package(ALPS CONFIG REQUIRED COMPONENTS applications solvers)`. The solver API is in `<alps/solvers.hpp>`. Python extensions that share ALPS objects with pyalps use its separate [downstream CMake package](python/pyalps/README.md#downstream-native-extensions).

Installation follows `GNUInstallDirs`. Unix SDKs continue to need their external Boost, HDF5 and numerical libraries after relocation. Redistributable Python wheels need auditwheel/delocate repair, as described in the [Python runtime guide](python/pyalps/README.md#native-runtime-layout).

### XML resources and tools

Installed programs use `share/alps/xml`. CTest supplies the source resources automatically; when running an uninstalled native program directly, set `ALPS_XML_PATH` to the absolute path of `src/alps/resources/`.

Unix application builds install `alps-xml`, which requires Python 3 and `xsltproc` on `PATH`. With the SDK's `bin` directory on `PATH`, use `alps-xml --help`, or, with your own result files:

```sh
alps-xml plot text results.plot.xml
alps-xml convert html simulation.out.xml --output results.html
alps-xml extract text plot-definition.xml task*.out.xml --output measurements.txt
```

Plot/extraction formats are `text`, `html`, `gnuplot`, `matplotlib` and `grace`; conversion supports `text` and `html`. The `xml` installation component includes the command and resources.

## CI coverage

Pull requests and master pushes report an aggregate `CI` check. The [PR workflow](.github/workflows/ci.yml) selects affected areas: static/helper checks, Linux native and editable Python tests, MPI integration, macOS SDK consumers, native sanitizers, and Linux wheel packaging. It uses binary dependency packages where available and compiler caches to reduce repeated builds.

Each area fingerprints its source, tests, build configuration and runner image. A successful matching fingerprint can skip a repeated check; missing history, unknown paths or uncertain parsing cause checks to run. Citation metadata and generated snapshots invalidate all areas. Use the workflow's manual `force` input to run every area regardless of recorded passes. See [change detection](ci/README.md) for details.

The [Compatibility workflow](.github/workflows/compatibility.yml) runs weekly and on demand. It checks minimum supported dependencies, extensive C++20/OpenMP tests, C++23, MPI-enabled Python, additional macOS platforms, installed tutorials, and the full manylinux/musllinux/macOS wheel matrix on CPython 3.11–3.14. PR packaging builds manylinux wheels and tests Python 3.11, 3.12 and 3.14; SDK consumers compile once per wheel platform, while newer interpreters check the reused abi3 artifact. Fresh-runner wheel imports and two-rank MPI adapter tests remain separate checks.

Release tags run Compatibility before publication through the [release workflow](.github/workflows/build_wheels.yml). Its existing filename is retained for the PyPI trusted publisher. Manual release validation does not publish. Local editable tests do not replace repaired-wheel or platform validation.

## Preparing a release

For maintainers:

1. Review the Unreleased changelog and migration notes. Give the release section its version and date, then add a new Unreleased section.
2. Update `cmake/ALPS_VERSION.txt`, shared by the SDK and Python package, to the numeric `X.Y.Z` release version. A final tag is `vX.Y.Z`; prereleases use `vX.Y.Z-alpha.N`, `-beta.N`, `-rc.N` or `-dev.N` while the file remains `X.Y.Z`. See [versioning](python/pyalps/README.md#versioning).
3. Validate the intended tag locally with Python ≥ 3.11, replacing `vX.Y.Z` below with the actual tag:

   ```sh
   python -m pip install packaging
   python .github/scripts/check_release_version.py --ref refs/tags/vX.Y.Z
   ```

4. Merge and validate the release commit before tagging. Tag-triggered workflows validate source builds, wheels and the source distribution before the PyPI upload job. Publication also requires the repository's `pypi` environment and trusted-publisher configuration; a tag in an arbitrary fork is not sufficient. Manual validation runs do not publish.

Keep published tags fixed. If a published release is wrong, correct the source and prepare a new version; do not use `skip-existing` to conceal mismatched artifacts.
