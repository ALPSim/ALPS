# Changelog

Release notes and migration guidance for ALPS. Changes awaiting release are collected under **Unreleased**; released versions receive a version number and date. The existing [looper history](src/apps/qmc/looper/ChangeLog) remains with that component.

## Unreleased

### Added

- Native Windows x64 and ARM64 builds with CMake presets, pinned vcpkg dependencies, and SDK installation containing the required runtime DLLs.
- A relocatable CMake SDK with `ALPS::alps`, `ALPS::headers`, and `ALPS::fortran`. Builds with applications also export executable targets and the solver libraries `ALPS::maxent`, `ALPS::cthyb`, and `ALPS::ctint`.
- The Unix `alps-xml` command for plot rendering, result conversion, and data extraction, including Python 3 Matplotlib output.
- A [tutorial guide](tutorials/README.md) that brings simulation tutorials, notebooks, and focused library examples together. Installed C++ and Fortran examples can be built against the SDK.

### Changed

- Python requires GIL-enabled CPython 3.12 or newer. Native wheels use the CPython 3.12 stable ABI (`cp312-abi3`), with one build per platform/architecture tested across Python 3.12–3.14. Downstream nanobind extensions exchanging ALPS objects must also enable `STABLE_ABI` and be rebuilt. Free-threaded Python is unsupported.
- Source builds require CMake 4.3 or newer and C++17. Dependencies and compiler requirements propagate through exported CMake targets.
- Source CI uses four routine configurations and fourteen weekly, release-tag, or manually requested configurations, including sanitizer coverage. Windows and stable-ABI wheel checks remain separate; workflow linting, verified Boost downloads, and test summaries make failures easier to diagnose.
- MPI is opt-in. Embedded builds default to the library alone. `BUILD_TESTING` controls the tests; `ALPS_BUILD_APPLICATIONS` controls applications and CLI tools; `ALPS_BUILD_EXTENSIVE_TESTS` adds the expensive graph and HDF5 tests. C++ examples are opt-in, and tutorials are a separate installation component.
- BLAS and LAPACK are required and use one LP64 ABI: 32-bit integers and lowercase symbols with a trailing underscore.
- Python bindings build as a separate `scikit-build-core` project against an installed C++ SDK. Python extensions reuse the SDK's solver libraries.
- Python builds require a shared SDK with the same numeric release version. `PYALPS_BUILD_SOLVERS` controls solver bindings; `PYALPS_BUNDLE_APPLICATIONS` independently controls bundled command-line programs. Required runtime libraries are always included, and extensions share a package-specific nanobind library.
- Installation follows `GNUInstallDirs`, supports relocated SDKs, and exposes XML resources through `ALPS_DATA_DIR` and the `ALPS_XML_PATH` environment override.
- The bundled XDR sources live under `third_party/xdr/`, with their license included in SDK and Python wheel installations. Public XDR header paths remain unchanged.
- Looper and hybridization solver references live under `tutorials/`, with current build instructions. The Looper reference uses Markdown; the hybridization manual retains its LaTeX source and bibliography.
- Repository paths now reflect their purpose:

  | Previous path | Current home |
  | --- | --- |
  | `applications/` | `src/apps/` |
  | `tool/` | CLI utilities in `src/tools/`; MaxEnt in `src/apps/maxent/` |
  | `lib/xml/` | `src/alps/resources/` |
  | `test/` | `tests/` |
  | `example/` | `tutorials/examples/` |
  | `bindings/boost/` | `third_party/boost_numeric_bindings/` |
  | `bindings/python/pyalps/` | `python/pyalps/`; native sources consolidated in `cpp/` |

  Retained CI and release helpers live in `.github/scripts/`.

### Removed

- Legacy CMake integration through `ALPS_USE_FILE`, `ALPS_LIBRARIES`, and dependency-variable aliases. Consumers link to the exported targets instead.
- The `ALPS_XML_PATH` CMake cache option and `alpsvars` environment scripts. The runtime environment override remains supported.
- The generated `pyalps_config.py` and its build-machine fallback paths. Python installations locate their own resources; external programs can be selected through `PATH` or `ALPS_BIN_PATH`.
- The individual XML shell commands, replaced by `alps-xml` subcommands.
- The SQLite archive application, its two XML archive converter frontends, standalone Alea analysis programs, and obsolete build and release scripts.
- The `src/boost/` compatibility directory, including obsolete codecvt templates, the Boost 1.61/1.62 Intel compiler patch, and an unused exception handler. Unreferenced Numeric Bindings headers and container adapters have also been removed.
- Redundant Boost.Accumulators copies under `src/alps/boost/` and the Mersenne Twister forwarding header. Include `<boost/random/mersenne_twister.hpp>` directly.
- Obsolete application-test scaffolding, the unused SSE2 tutorial, a committed editor cache, and the generated hybridization PDF. The active XML CLI fixture now lives beside its test.

### Fixed

- Default observable data initializes its thermalization state. Integer ranges retain valid bounds across the entire integer domain; `size()` returns `std::uintmax_t` and throws `std::overflow_error` when the count cannot fit that type.
- Simulation-help issues route Worm and SSE to their respective maintainers without an external form parser. DMFT startup notices preserve scientific credit while removing duplicate citation requests and mandatory-citation wording.
- Windows ARM64 numerical builds use a compatible BLAS/LAPACK provider instead of mixing incompatible CLAPACK and OpenBLAS return conventions.
- ALPS uses the upstream Spirit Classic headers and namespace directly, so including Boost Spirit before ALPS no longer depends on an ALPS compatibility macro.
- Installed examples include their parameter files and fixed HDF5 input; their tests locate SDK resources and Windows runtime DLLs.
- Editable Python installations use live Python sources while locating compiled extensions and runtime resources in their installed directory. A bindings-only Windows package no longer masks separately installed command-line programs with its DLL directory.

### Migration

1. Use a fresh build directory and installation prefix when upgrading from an older layout. `cmake --install` does not remove obsolete headers or executables.
2. Replace legacy CMake variables and use files with `find_package(ALPS CONFIG REQUIRED)` and `target_link_libraries(... ALPS::alps)`. See [SDK usage](CONTRIBUTING.md#consuming-the-c-sdk) for additional targets.
3. Build Python bindings against the installed SDK. Replace the Python build option `PYALPS_BUILD_APPLICATIONS` with `PYALPS_BUILD_SOLVERS`. Downstream native Python extensions use the [CMake package supplied by pyalps](python/pyalps/README.md#downstream-native-extensions).
4. Update XML-tool invocations to the [`alps-xml` subcommands](CONTRIBUTING.md#xml-command-line-tools). Run Python tutorials with an interpreter from an environment containing pyalps.
5. Update source paths in local scripts using the table above. Current build switches and platform setup are documented in [the contributor guide](CONTRIBUTING.md#getting-started-with-the-code).
