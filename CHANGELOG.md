# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Changed

- Require CMake 3.27 or newer and an externally installed Boost 1.76 or newer with CMake packages. The SDK requires C++17/C11 compilers, HDF5's C library, and LP64 BLAS/LAPACK; bundled Boost builds and alternate numerical integer/symbol ABIs are no longer supported.
- Export CMake targets for the installed SDK, applications and solver libraries. Downstream C++ projects link `ALPS::alps`; Python extensions sharing pyalps objects use `pyalps::runtime`. MaxEnt, CT-HYB and CT-INT Python wrappers link the SDK's solver libraries instead of compiling their implementations again.
- Make MPI opt-in with `ALPS_ENABLE_MPI=ON`. Standalone builds enable applications and native tests by default; embedded `add_subdirectory` builds default to the library alone. The default SDK uses shared libraries, as required by the Python bindings.
- Support editable Python development against an installed SDK, with packaged native runtime resources and explicit downstream ABI checks. Redistributable wheels bundle native dependencies through platform repair and are tested after installation on fresh runners.
- Organize tutorials into numbered topic directories, with introductory scripts directly in `tutorials/01-intro/`. `tutorials/00-examples/` is an independent API reference collection. Tutorial sources are installed only when the `tutorials` component is requested. See the [tutorial index](tutorials/README.md).

### Removed and migration

The old build interfaces are removed without compatibility aliases. Reconfigure existing build scripts and downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

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
| XML definitions and stylesheets in `lib/xml/` | `src/alps/resources/` in the source tree; installed under `share/alps/xml/` |
| `tutorials/alpsize-*`, `tutorials/code-*`, `tutorials/ngs/` | `tutorials/08-alpsize/*`, `tutorials/09-code/*`, `tutorials/10-ngs/`; the Python export example is now in `python/pyalps/examples/ising/` |
| Root `ALPS_VERSION.txt` | `cmake/ALPS_VERSION.txt`, still shared by the SDK and Python package |

### Fixed

- Generate XML plotting scripts compatible with Python 3 through the `alps-xml` CLI.
- Correct the Heisenberg tutorial's vector dot product and allow its vector implementation to compile without x86 SIMD support.
- Remove undefined behavior in XDR callbacks and integer-range boundary checks exposed by sanitizers.
