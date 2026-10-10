# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Fixed

- Release HDF5 scalar string buffers after reads, including failed conversions, and owned observables when clearing an `ObservableSet`.
- Correct scalar/vector Monte Carlo result negation and keep uncertainties nonnegative under negative scaling or reciprocal arithmetic. Preserve empty histogram ranges, avoid invalid access during empty vector/valarray conversions, and create private temporary files on Unix independently of the caller's umask.

### Changed

- Remove obsolete build, dependency-bootstrap and installer scripts. Replace historical XML shell wrappers with `alps-xml plot`, `alps-xml convert` and `alps-xml extract`, and update notebook commands. Keep installed `txt2archive`, `xml2archive` and the optional SQLite archive tool under `src/tools/archive`.

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
- Organize tutorials into numbered topic directories, with introductory scripts directly in `tutorials/01-intro/`. `tutorials/examples/` is an independent API reference collection. Tutorial sources are installed only when the `tutorials` component is requested. See the [tutorial index](tutorials/README.md).

### Migration

The existing native test framework, public and vendored headers, and differential audit tooling are retained. Tutorial consumers use the exported SDK targets; the obsolete Makefile and CMake consumer interfaces are removed. Configure downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

Numerical matrix/vector persistence now requires an explicit adapter: include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` and link `ALPS::numeric_io`. The matrix umbrella `<alps/numeric/matrix.hpp>` no longer includes HDF5 automatically. Diagonal matrices and deprecated BLAS matrix/vector classes retain their archive `save`/`load` members without requiring HDF5 headers in the numerical interfaces. Numerical matrix `write_xml` and XML insertion remain available through `<alps/xml/matrix.hpp>` and `ALPS::numeric_xml`. Header ownership moves preserve public include spellings; consumers of persistence and XML output must include the corresponding adapters.

| Previous interface or location | Replacement |
| --- | --- |
| `wheel-deps` preset | `distribution` preset; update `ALPS_DIR` to the new installation |
| `ALPS_BUILD_TESTS` | `ALPS_BUILD_TESTING`; a parent's `BUILD_TESTING` does not control ALPS tests |
| `ALPS_BUILD_LIBS_ONLY=ON` | `ALPS_BUILD_APPLICATIONS=OFF`; the `sdk` preset also disables tests |
| `ALPS_BUILD_EXAMPLES` | Build the standalone projects in `tutorials/examples/` against an installed SDK |
| `ALPS_BUILD_FORTRAN` | The SDK always includes `ALPS::fortran`; build Fortran examples separately with a Fortran compiler |
| `ALPS_INSTALL_HEADERS` | SDK headers are always installed |
| `ALPS_INCLUDE_TUTORIALS` | `cmake --install <build-dir> --component tutorials` |
| `ALPS_ENABLE_OPENMP_WORKER` | `ALPS_ENABLE_OPENMP`; choose the worker thread count at runtime |
| `PYALPS_BUILD_APPLICATIONS` | `PYALPS_BUILD_SOLVERS`; `PYALPS_BUNDLE_APPLICATIONS` still controls bundled executables separately |
| `ALPS_USE_SYSTEM_BOOST`, `Boost_ROOT_DIR` / bundled Boost discovery | External Boost packages are always required; locate them through `Boost_ROOT` or `CMAKE_PREFIX_PATH` |
| `LAPACK_64_BIT`, alternate `BIND_FORTRAN_*` ABIs | Use LP64 BLAS/LAPACK with lowercase, trailing-underscore symbols |
| `UseALPS.cmake`, `include.mk`, `alpsvars` scripts | Imported SDK targets and explicit installation paths; add the installed `bin` directory to `PATH` |
| `UsePyALPS.cmake`, `<alps/ngs/detail/export_sim_to_python.hpp>` | [pyalps downstream CMake package](python/pyalps/README.md#downstream-native-extensions) and `<pyalps/export_simulation.hpp>` |
| Historical XML shell tools | Retained alongside `alps-xml`; see [XML tools](CONTRIBUTING.md#xml-resources-and-tools) |
| Top-level `import mpi` compatibility module | `import pyalps.mpi`; install the `mpi` extra for mpi4py |
| `applications/`, `tool/`, `test/`, `example/` | `src/apps/`, `src/tools/`, `tests/`, `tutorials/examples/` |
| Loose subsystem trees and transitional `src/alps/{common,runtime}/` | Semantic `src/alps/<module>/{include,src,tests}` ownership; see the [module map](src/alps/README.md#source-ownership) |
| `src/ietl/` | `src/alps/ietl/include/ietl/`; public `<ietl/...>` includes are unchanged |
| `src/alps/config.h.in`, `src/alps/version.h.in`, transitional `src/alps/common/config/` | `cmake/config/` configuration templates |
| Build-tree `src/alps/` generated headers | `<build-dir>/generated/include/alps/`; consumers use exported CMake targets |
| XML definitions and stylesheets in `lib/xml/` | `src/alps/resources/` in the source tree; installed under `share/alps/xml/` |
| `tutorials/alpsize-*`, `tutorials/code-*`, `tutorials/ngs/` | `tutorials/08-alpsize/*`, `tutorials/09-code/*`, `tutorials/10-ngs/`; the Python export example is now in `python/pyalps/examples/ising/` |
| Root `ALPS_VERSION.txt` | `cmake/ALPS_VERSION.txt`, still shared by the SDK and Python package |

### Fixed

- Initialize Monte Carlo completion scheduling before first use.
- Read integer checkpoint parameters at their stored width and reject out-of-range conversion without changing existing values.
- Preserve `txt2archive`, `xml2archive`, and the optional SQLite `archive` program; correct malformed-input handling.

- Report MaxEnt CLI help and input errors with normal exit codes instead of continuing into an empty input or aborting on an exception. Valid scientific runs are unchanged.
- Reject unsupported native C++ parameter checkpoint datatypes with the dataset path in the diagnostic, instead of silently substituting zero. Existing supported types and custom readers retain their decoding behavior; a failed parameter reload preserves the previous values.
- Generate XML plotting scripts compatible with Python 3 through the `alps-xml` CLI.
- Correct the Heisenberg tutorial's vector dot product and allow its vector implementation to compile without x86 SIMD support.
- Remove undefined behavior in XDR callbacks and integer-range boundary checks exposed by sanitizers.
