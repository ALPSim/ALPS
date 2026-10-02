# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Changed

- Export `ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` as interface targets. Foundation libraries and interfaces own their header sets without inheriting aggregate `ALPS::headers`. Separate array storage from mathematical helpers and numerical types from HDF5 adapters, removing the foundation include cycle; the expression/older-parameters cycle remains. Numerical algorithms and the moved headers' public include names are preserved, with no ALPSCore implementation imported. Params conversion headers move to `params/adapters/include/` and remain exposed through the aggregate interface with implementations in `ALPS::alps`.
- Group `src/tools/` commands by parameter preparation, lattice export, scheduler formats, Parapack, diagnostics and XML transformations. Keep historical inactive tools under explicit owners without enabling them. Executable names, installation components and source contents are preserved; `pconfig` now links only the utilities component.
- Organize `src/alps/` by semantic ownership, including numerical/container helpers, XML, older parameters, models, observables and execution, replacing the transitional `common/` and `runtime/` groups. Modules own their headers, sources and local tests; configuration templates live in `cmake/config/`. Explicit file sets preserve public include names, and generated headers use `<build-dir>/generated/include/alps/`. A CMake-generated module manifest supports ownership and include-dependency checks. Remaining simulation modules use the aggregate header interface; physical owners are not all independent libraries.
- Export `ALPS::xml` for XML parsing/output and `ALPS::cli` for the existing `mcoptions` and `parseargs` command-line grammars. MaxEnt's executable now links `ALPS::cli` instead of `ALPS::alps`; its HDF5 input, command-line behavior and scientific calculations are unchanged. Move `<alps/plot.h>` to the `plotting` header module because it also uses older parameters. Public include names remain stable, and `ALPS::headers` remains the aggregate compile interface. Python runtime packaging includes both libraries; rebuild downstream binaries after the split.
- Export `ALPS::osiris` for dump/XDR serialization and process/communication state. MaxEnt's solver links the foundations, Osiris and numerical providers without the simulation runtime, preserving its numerical calculations and stop-callback behavior. Python runtime packaging includes Osiris. Rebuild downstream binaries after the library split.
- Export `ALPS::params` as a separate typed parameter library with component-owned exports. Isolate the existing parameter-file constructor, XML input and older `Parameters` conversion in adapters owned by `ALPS::alps`; their source APIs and parsing behavior are preserved. Python packages carry the shared runtime components. Rebuild downstream binaries after the split.
- Export `ALPS::hdf5` as a separate archive library with its own symbol exports and NGS signal cleanup. It links utilities and HDF5/Boost dependencies without the simulation runtime; `ALPS::alps` links it transitively. Python packages include one copy of this runtime component. Rebuild downstream binaries after the split; archive formats and public include paths are unchanged.
- Export `ALPS::utilities` as a separate library with its own symbol exports. `ALPS::alps` links it transitively, and Python packages carry this runtime component. Rebuild downstream binaries after this library split; source include paths and utility APIs are unchanged.
- Group utilities, HDF5 and NGS parameter headers, implementations and tests by module under `src/alps/`; separate MaxEnt's implementation, CLI and tests under `src/apps/maxent/`. Public include paths and exported library targets are unchanged. See the [module layout](src/alps/README.md) for the ALPSCore reconciliation boundaries.

- Require CMake 3.27 or newer and an externally installed Boost 1.76 or newer with CMake packages. The SDK requires C++17/C11 compilers, HDF5's C library, and LP64 BLAS/LAPACK; bundled Boost builds and alternate numerical integer/symbol ABIs are no longer supported.
- Export CMake targets for the installed SDK, applications and solver libraries. Downstream C++ projects link `ALPS::alps`; Python extensions sharing pyalps objects use `pyalps::runtime`. MaxEnt, CT-HYB and CT-INT Python wrappers link the SDK's solver libraries instead of compiling their implementations again.
- Make MPI opt-in with `ALPS_ENABLE_MPI=ON`. Standalone builds enable applications and native tests by default; embedded `add_subdirectory` builds default to the library alone. The default SDK uses shared libraries, as required by the Python bindings.
- Support editable Python development against an installed SDK, with packaged native runtime resources and explicit downstream ABI checks. Redistributable wheels bundle native dependencies through platform repair and are tested after installation on fresh runners.
- Organize tutorials into numbered topic directories, with introductory scripts directly in `tutorials/01-intro/`. `tutorials/00-examples/` is an independent API reference collection. Tutorial sources are installed only when the `tutorials` component is requested. See the [tutorial index](tutorials/README.md).

### Removed and migration

The old build interfaces are removed without compatibility aliases. Reconfigure existing build scripts and downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

Numerical persistence now requires an explicit adapter: include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` and link `ALPS::numeric_io`. The matrix umbrella `<alps/numeric/matrix.hpp>` no longer includes HDF5. HDF5 `save`/`load` members on diagonal matrices and deprecated BLAS matrix/vector classes, plus numerical matrix `write_xml` and XML insertion overloads, were unused in this repository and are removed without replacement adapters. Header ownership moves preserve public include spellings, but these method removals and explicit adapter requirements are source API changes.

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
