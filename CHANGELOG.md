# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Changed

- Fix defects exposed by sanitizer coverage: reclaim scalar HDF5 variable-length
  strings, release owned observables when clearing an observable set, and handle
  empty vector/valarray conversions without indexing nonexistent elements.
- Standardize native runtime tests on GoogleTest with individually discoverable
  CTest cases, isolated fixtures, numerical assertions, and preserved historical
  serialization contracts. Add development, MPI, extensive, and sanitizer test
  presets; see [the testing guide](tests/README.md). GoogleTest is required only
  when building tests and is never installed with the SDK.
- Run standalone NumPy tutorial checks independently of native builds. Use
  representative PR configurations and the full supported matrix for shared
  build changes, scheduled validation and releases; retain installed-artifact,
  multi-rank MPI and sanitizer checks with machine-readable test reports.

- Require HDF5 1.10.5 or newer for source builds and installed SDK consumers, retaining compatibility with the system package used by the manylinux_2_28 wheel build.
- Require Python 3.11 or newer for pyalps. Keep a separate CPython 3.11 wheel and use one `cp312-abi3` wheel for Python 3.12 and newer. Downstream nanobind extensions must pass `STABLE_ABI` to share pyalps types on Python 3.12+.
- Export `ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` as interface targets with their own header sets and dependencies. Separate container storage, numerical algorithms and HDF5 adapters; foundations no longer inherit aggregate `ALPS::headers`. Numerical algorithms and public include names are preserved.
- Organize `src/alps/` by responsibility, with module-local headers, sources and tests. Configuration templates live in `cmake/config/`, and generated headers use `<build-dir>/generated/include/alps/`. Explicit header file sets preserve public include names. See the [module layout](src/alps/README.md).
- Export independently linkable runtime components `ALPS::utilities`, `ALPS::hdf5`, `ALPS::params`, `ALPS::osiris`, `ALPS::xml` and `ALPS::cli`. `ALPS::alps` links them transitively; Python packages carry one copy of each component. Params text/XML and older `Parameters` conversion adapters remain in `ALPS::alps`, with unchanged public header names. Archive formats and existing parsing behavior are preserved. Rebuild downstream binaries after the library splits.
- MaxEnt's solver and executable link the foundation components without `ALPS::alps`. The executable uses `ALPS::cli` and reads typed params directly from HDF5. Scientific calculations and stop-callback behavior are preserved.
- Group `src/tools/` commands by responsibility. Keep historical inactive sources without enabling them; executable names and installation components are preserved. `pconfig` now links only utilities.
- Require CMake 3.27 or newer and an externally installed Boost 1.76 or newer with CMake packages. The SDK requires C++17/C11 compilers, HDF5's C library, and LP64 BLAS/LAPACK; bundled Boost builds and alternate numerical integer/symbol ABIs are no longer supported.
- Export CMake targets for the installed SDK, applications and solver libraries. Downstream C++ projects link `ALPS::alps`; Python extensions sharing pyalps objects use `pyalps::runtime`. MaxEnt, CT-HYB and CT-INT Python wrappers link the SDK's solver libraries instead of compiling their implementations again.
- Make MPI opt-in with `ALPS_ENABLE_MPI=ON`. Standalone builds enable applications and native tests by default; embedded `add_subdirectory` builds default to the library alone. The default SDK uses shared libraries, as required by the Python bindings.
- Support editable Python development against an installed SDK, with packaged native runtime resources and explicit downstream ABI checks. Redistributable wheels bundle native dependencies through platform repair and are tested after installation on fresh runners.
- Organize tutorials into numbered topic directories, with introductory scripts directly in `tutorials/01-intro/`. `tutorials/00-examples/` is an independent API reference collection. Tutorial sources are installed only when the `tutorials` component is requested. See the [tutorial index](tutorials/README.md).

### Removed and migration

The old build interfaces are removed without compatibility aliases. Reconfigure existing build scripts and downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

Remove the unused, incomplete headers under `<alps/ngs/accumulator/deprecated/>` and their unbuilt sources. The active `<alps/ngs/accumulator.hpp>` API remains available. Also remove the disconnected prototype accumulator tests, orphan Boost inspection tools, obsolete line-count script and unused IETL configuration template. The historical `<alps/ngs/scheduler/proto/tcpserver.hpp>` include now forwards to the canonical scheduler header.

Remove the inactive SQLite result archive, uninstalled Alea mean/variance frontends and historical XML shell wrappers superseded by `alps-xml`. Retire the one-off Boost.Python comparison harness; historical checkpoint loading is covered by a frozen fixture in the regular Python suite.

Remove the obsolete numerical containers under `<alps/numeric/deprecated/>` and `<alps/numeric/detail/deprecated/>`, and `<alps/numeric/resizeable_vector.h>`. Use the active `<alps/numeric/matrix.hpp>` and `<alps/numeric/matrix/vector.hpp>` interfaces instead. The hybridization solver's separate matrix implementation is unchanged. Prune unused vendored Numeric Bindings LAPACK wrappers after replacing ALPS's umbrella includes with specific dependencies; downstream code using the removed vendor headers must supply its own Numeric Bindings installation or migrate to ALPS's numerical interfaces.

Numerical matrix/vector persistence now requires an explicit adapter: include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` and link `ALPS::numeric_io`. The matrix umbrella `<alps/numeric/matrix.hpp>` no longer includes HDF5 automatically. Diagonal matrices retain their archive `save`/`load` members without requiring HDF5 headers in the numerical interfaces. Numerical matrix `write_xml` and XML insertion overloads are removed in 3.0 to keep numerical interfaces independent of XML; callers must handle XML output explicitly. Header ownership moves preserve public include spellings, but the explicit adapter requirement and XML method removals are source API changes.

| Previous interface or location | Replacement |
| --- | --- |
| `wheel-deps` preset | `distribution` preset; update `ALPS_DIR` to the new installation |
| `ALPS_BUILD_TESTS` | `ALPS_BUILD_TESTING`; a parent's `BUILD_TESTING` does not control ALPS tests |
| `ALPS_BUILD_LIBS_ONLY=ON` | `ALPS_BUILD_APPLICATIONS=OFF`; the `sdk` preset also disables tests |
| `ALPS_BUILD_EXAMPLES` | Build the standalone projects in `tutorials/00-examples/` against an installed SDK |
| `ALPS_BUILD_FORTRAN` | The SDK always includes `ALPS::fortran`; build Fortran examples separately with a Fortran compiler |
| `ALPS_INSTALL_HEADERS` | SDK headers are always installed |
| `ALPS_INCLUDE_TUTORIALS` | `cmake --install <build-dir> --component tutorials` |
| `ALPS_ENABLE_OPENMP_WORKER` | `ALPS_ENABLE_OPENMP`; choose the worker thread count at runtime |
| `PYALPS_BUILD_APPLICATIONS` | `PYALPS_BUILD_SOLVERS`; `PYALPS_BUNDLE_APPLICATIONS` still controls bundled executables separately |
| `ALPS_USE_SYSTEM_BOOST`, `Boost_ROOT_DIR` / bundled Boost discovery | External Boost packages are always required; locate them through `Boost_ROOT` or `CMAKE_PREFIX_PATH` |
| `LAPACK_64_BIT`, alternate `BIND_FORTRAN_*` ABIs | Use LP64 BLAS/LAPACK with lowercase, trailing-underscore symbols |
| `UseALPS.cmake`, `include.mk`, `alpsvars` scripts | Imported SDK targets and explicit installation paths; add the installed `bin` directory to `PATH` |
| `UsePyALPS.cmake`, `<alps/ngs/detail/export_sim_to_python.hpp>` | [pyalps downstream CMake package](python/pyalps/README.md#downstream-native-extensions) and `<pyalps/export_simulation.hpp>` |
| `plot2*`, `convert2html`, `convert2text`, `extract*` shell tools | `alps-xml plot`, `alps-xml convert`, `alps-xml extract`; see [XML tools](CONTRIBUTING.md#xml-resources-and-tools) for formats and dependencies |
| Top-level `import mpi` compatibility module | `import pyalps.mpi`; install the `mpi` extra for mpi4py |
| `applications/`, `tool/`, `test/`, `example/` | `src/apps/`, `src/tools/`, `tests/`, `tutorials/00-examples/` |
| Loose subsystem trees and transitional `src/alps/{common,runtime}/` | Semantic `src/alps/<module>/{include,src,tests}` ownership; see the [module map](src/alps/README.md#source-ownership) |
| `src/ietl/` | `src/alps/ietl/include/ietl/`; public `<ietl/...>` includes are unchanged |
| `src/alps/config.h.in`, `src/alps/version.h.in`, transitional `src/alps/common/config/` | `cmake/config/` configuration templates |
| Build-tree `src/alps/` generated headers | `<build-dir>/generated/include/alps/`; consumers use exported CMake targets |
| XML definitions and stylesheets in `lib/xml/` | `src/alps/resources/` in the source tree; installed under `share/alps/xml/` |
| `tutorials/alpsize-*`, `tutorials/code-*`, `tutorials/ngs/` | `tutorials/08-alpsize/*`, `tutorials/09-code/*`, `tutorials/10-ngs/`; the Python export example is now in `python/pyalps/examples/ising/` |
| Root `ALPS_VERSION.txt` | `cmake/ALPS_VERSION.txt`, still shared by the SDK and Python package |

### Fixed

- Report MaxEnt CLI help and input errors with normal exit codes instead of continuing into an empty input or aborting on an exception. Valid scientific runs are unchanged.
- Reject unsupported native C++ parameter checkpoint datatypes with the dataset path in the diagnostic, instead of silently substituting zero. Existing supported types and custom readers retain their decoding behavior; a failed parameter reload preserves the previous values.
- Generate XML plotting scripts compatible with Python 3 through the `alps-xml` CLI.
- Correct the Heisenberg tutorial's vector dot product and allow its vector implementation to compile without x86 SIMD support.
- Remove undefined behavior in XDR callbacks and integer-range boundary checks exposed by sanitizers.
