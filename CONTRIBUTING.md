# Contributing to ALPS

Thank you for your interest in ALPS (Algorithms and Libraries for Physics Simulations).
ALPS is a community-driven, open-source ecosystem for numerical simulations of correlated quantum systems.
Contributions at every level — from a one-line bug report to a new simulation method — are welcome and valued.

## Table of contents

- [Ways to contribute](#ways-to-contribute)
- [Reporting bugs and requesting features](#reporting-bugs-and-requesting-features)
- [Getting started with the code](#getting-started-with-the-code)
- [Making a change](#making-a-change)
- [Provenance and scientific credit](#provenance-and-scientific-credit)
- [Submitting a pull request](#submitting-a-pull-request)
- [Build reference](#build-reference)
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
For third-party material, also follow the [provenance guidance](#provenance-and-scientific-credit) below.

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

- CMake ≥ 3.27, Ninja for the bundled presets, and C++17/C11 compilers such as GCC or Clang.
- Boost ≥ 1.76 with its compiled libraries and CMake packages, HDF5 ≥ 1.10.5 (C library), and LP64 BLAS/LAPACK. Use serial HDF5 for the default MPI-disabled build; see [numerical libraries](#numerical-libraries).
- For Python development: GIL-enabled CPython ≥ 3.10 in a writable Python environment. Pip installs NumPy, SciPy and Matplotlib with pyalps. Free-threaded Python is unsupported.
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

CMake and Ninja can also come from system packages or [official CMake binaries](https://cmake.org/download/). A generator other than Ninja can be selected with a plain CMake invocation instead of a preset.

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

The Python bindings are a separate `scikit-build-core` project that builds
against an installed ALPS C++ SDK; see the
[`pyalps` build instructions](python/pyalps/README.md).

### Citation maintenance

Citation metadata and application rules live in `CITATION.cff` and `CITATIONS.yaml`. See [citation maintenance](.github/scripts/citations/README.md) for generation and validation. Native builds use checked-in generated citation data and do not require Python.

Editing that data requires Python ≥ 3.9 with `PyYAML` and `jsonschema`: `python -m pip install -r .github/scripts/citations/requirements.txt`, then `python .github/scripts/generate_citations.py --regenerate`. Commit the generated files alongside the authorities; CI checks that they agree. For the optional citation tests, select an interpreter with `-DALPS_CITATION_PYTHON=/path/to/python`.

### Run the tests

From the repository root:
```sh
ctest --preset default
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
| `BUILD_SHARED_LIBS` | `ON` when unset | Build shared SDK libraries |
| `ALPS_ENABLE_MPI` | `OFF` | Enable MPI; requires matching MPI and Boost.MPI installations |
| `ALPS_ENABLE_OPENMP` | `OFF` | Enable OpenMP, including worker scheduling |
| `ALPS_BUILD_EXTENSIVE_TESTS` | `OFF` | Add expensive graph and HDF5 type-matrix tests when testing is enabled |

For example, configure with `cmake --preset default -DALPS_ENABLE_OPENMP=ON`. The `sdk` preset disables applications and tests; `distribution` disables tests and its build preset installs automatically. When embedding ALPS with `add_subdirectory`, applications and tests default to `OFF`. MPI remains opt-in. Headers and the C++ Fortran bridge are always part of the SDK; building that bridge needs no Fortran compiler.

Examples build separately against an installed SDK; see the [C++ and Fortran example instructions](tutorials/00-examples/README.md). Fortran tutorials that call OpenMP also need a Fortran OpenMP runtime. To install tutorial sources under `share/alps/tutorials`, run `cmake --install _build/default --component tutorials` after installing the SDK.

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

The imported target carries headers, C++17 requirements, compile definitions and transitive dependencies. Use a compiler and configuration compatible with the SDK's ABI. `ALPS::headers` exposes the compile interface without linking; `ALPS::fortran` supplies the C++ Fortran bridge and its GNU Fortran compatibility flag. Programs using only utilities can request `find_package(ALPS CONFIG REQUIRED COMPONENTS utilities)` and link `ALPS::utilities`; MPI-enabled SDKs also carry the runtime dependencies of `<alps/ngs/boost_mpi.hpp>`. Archive-only programs can similarly request the `hdf5` component and link `ALPS::hdf5`, which also owns archive signal cleanup. Typed parameter programs can request the `params` component and link `ALPS::params`. XML parsing/output and command-line parsing are available through the `xml` and `cli` components and targets `ALPS::xml` and `ALPS::cli`. Legacy text loading uses `alps::params_from_file(path)` from `<alps/ngs/params_from_file.hpp>`; this free function and the text/XML conversion adapters require `ALPS::alps`. Package discovery still checks the SDK's complete dependency set; component-specific configuration is future work.

Use `ALPS::containers` for array storage and `ALPS::numerics` for matrix/vector algorithms and array mathematics, including the existing `<alps/multi_array.hpp>` umbrella. Numerical archive consumers link `ALPS::numeric_io` and explicitly include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>`. For example:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS numeric_io)
target_link_libraries(my_simulation PRIVATE ALPS::numeric_io)
```

`<alps/numeric/matrix.hpp>` no longer includes its HDF5 adapter automatically.
For matrix XML output, include `<alps/xml/matrix.hpp>` and link
`ALPS::numeric_xml` (component `numeric_xml`). Both `matrix.write_xml(xml)` and
`xml << matrix` retain the historical `MATRIX`/`ROW`/`ELEMENT` representation.
The numerical interfaces themselves do not depend on XML or HDF5.

An SDK built with applications also exports executable targets such as `ALPS::spinmc` and the solver libraries `ALPS::maxent`, `ALPS::cthyb` and `ALPS::ctint`. Require them with `find_package(ALPS CONFIG REQUIRED COMPONENTS applications solvers)`. The solver API is in `<alps/solvers.hpp>`.

Installation follows `GNUInstallDirs`. Unix SDKs continue to need their external Boost, HDF5 and numerical libraries after relocation.

The component reorganization and `<alps/solvers.hpp>` API are development changes
for a future release after 3.0. ALPSCore HDF5 and parameter reconciliation remains
separate work. `ALPS::params` currently
provides NGS params, not ALPSCore params; the solver interface uses that type.
Review downstream compatibility against the released 3.0. ALPS and ALPSCore
have overlapping `alps/` headers and C++ symbols: use separate prefixes and
processes while comparing implementations. A port must select one implementation
per interface, starting with HDF5 format compatibility and then parameter types.

### XML resources and tools

Installed programs use `share/alps/xml`. CTest supplies the source resources automatically; when running an uninstalled native program directly, set `ALPS_XML_PATH` to the absolute path of `src/alps/resources/`.

Unix application builds install `alps-xml`, which requires Python 3 and `xsltproc` on `PATH`. With the SDK's `bin` directory on `PATH`, use `alps-xml --help`, or, with your own result files:

```sh
alps-xml plot text results.plot.xml
alps-xml convert html simulation.out.xml --output results.html
alps-xml extract text plot-definition.xml task*.out.xml --output measurements.txt
```

Plot/extraction formats are `text`, `html`, `gnuplot`, `matplotlib` and `grace`; conversion supports `text` and `html`. The `xml` installation component includes the command and resources.

Application builds also retain the archive-producing commands in the `tools`
installation component:

```sh
txt2archive --xaxis T --yaxis Energy --with-error measurements.txt > results.archive.xml
xml2archive simulation.in.xml > results.archive.xml
```

`txt2archive` reads x/y columns, with an optional third error column;
`xml2archive` collects a master JOB XML file's task results. These operations
produce archive XML rather than render it. If SQLite development headers and
libraries are available at configure time, the `archive` tool is also built.
It retains its existing SQLite schema and `install`, `append`, `rebuild`, `list`
and `plot` commands; existing databases can be reopened without conversion.
Use `archive --help` for its options. The SDK libraries do not depend on SQLite.

## Preparing a release

Keep `cmake/ALPS_VERSION.txt` (the C++ SDK version) and the root
`ALPS_VERSION.txt` (the Python package version source) equal before creating
a release tag. Both files contain the numeric `X.Y.Z` core. A final release
uses the tag `vX.Y.Z`; prereleases use `vX.Y.Z-beta.N`, `-alpha.N`, `-rc.N`
or `-dev.N`. Python packaging derives the corresponding PEP 440 version.

Validate the intended tag locally using Python 3.11 or newer:

```bash
python -m pip install packaging
python .github/scripts/check_release_version.py --ref refs/tags/vX.Y.Z
```

The packaging workflow checks these versions before building and checks every
wheel and source distribution, including its embedded metadata, before upload.
Tag pushes publish the full release to PyPI, including CPython 3.10–3.14 wheels.
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

- CMake ≥ 3.27 features are acceptable.
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
