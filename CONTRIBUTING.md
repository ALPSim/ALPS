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

- CMake ≥ 3.22
- A C++17-capable compiler (GCC, Clang, Intel, or Fujitsu)
- Boost (downloaded automatically during configuration; or use a system install with `-DALPS_USE_SYSTEM_BOOST=ON`)
- For Fortran bindings: gfortran (or compatible Fortran compiler)
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
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

Alternatively, use the bundled CMake preset:
```bash
cmake --preset default
cmake --build --preset default
```

The Python bindings are a separate `scikit-build-core` project that builds
against an installed ALPS C++ SDK; see the
[`pyalps` build instructions](bindings/python/pyalps/README.md).

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
Only a tag push publishes to PyPI. Merge and validate the release commit before
tagging it.

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

- CMake ≥ 3.22 features are acceptable.
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
