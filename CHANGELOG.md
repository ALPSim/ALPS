# Changelog

User-facing changes and migration notes are recorded here, starting with the build modernization. Earlier releases are not yet catalogued in this file.

## Unreleased

### Fixed

- Store an empty `std::map` in an HDF5 archive as an empty group, so it reads back as an empty map instead of a missing path. Older readers load such files as empty maps too.
- Resolve HDF5 attribute paths written as `path@name` (ALPSCore's spelling) and attributes of the root group (`/@name`, or `@name` in the root context). These used to drop the character before `@` or fail; ALPS's `path/@name` spelling names the same object as before.
- Release HDF5 scalar string buffers after reads, including failed conversions, and owned observables when clearing an `ObservableSet`.
- Correct scalar/vector Monte Carlo result negation and keep uncertainties nonnegative under negative scaling or reciprocal arithmetic. Preserve empty histogram ranges, avoid invalid access during empty vector/valarray conversions, and create private temporary files on Unix independently of the caller's umask.
- Initialize Monte Carlo completion scheduling before first use.
- Read integer checkpoint parameters at their stored width and reject out-of-range conversion without changing existing values.
- Preserve `txt2archive`, `xml2archive`, and the optional SQLite `archive` program; correct malformed-input handling.
- Report MaxEnt CLI help and input errors with normal exit codes instead of continuing into an empty input or aborting on an exception. Valid scientific runs are unchanged.
- Reject unsupported native C++ parameter checkpoint datatypes with the dataset path in the diagnostic, instead of silently substituting zero. Existing supported types and custom readers retain their decoding behavior; a failed parameter reload preserves the previous values.
- Generate XML plotting scripts compatible with Python 3 through the `alps-xml` CLI.
- Correct the Heisenberg tutorial's vector dot product and allow its vector implementation to compile without x86 SIMD support.
- Remove undefined behavior in XDR callbacks and integer-range boundary checks exposed by sanitizers.
- Loading a `std::map` from HDF5 through an empty path inside an archive context, as a class's `load(archive&)` does, read the entries below the root instead of below the context. Entry paths are now completed against the context, as saving already did and as ALPSCore does.
- Saving a `std::map` at an HDF5 archive's root tried to delete the root group and failed. It now writes its entries beside the root's other children.

### Changed

- Consolidate native runtime tests under GoogleTest with individual CTest cases, historical numerical and serialization references, isolated temporary files, and explicit MPI/sanitizer registration. GoogleTest is a test-only dependency; builds with tests disabled do not acquire or install it.
- Keep Python tests and fixtures alongside pyalps and run standalone numerical tutorial tests independently. Retire the optional Boost.Python differential build audit after retaining application checks and focused regressions; see [the coverage assessment](tests/python-migration.md) for historical compatibility limits.

- Require HDF5 1.10.5 or newer for source builds and installed SDK consumers, retaining compatibility with the system package used by the manylinux_2_28 wheel build.
- Support Python 3.10 or newer for pyalps. Keep separate CPython 3.10 and 3.11 wheels and use one `cp312-abi3` wheel for Python 3.12 and newer. Downstream nanobind extensions must pass `STABLE_ABI` to share pyalps types on Python 3.12+.
- Add `mpi`, `openmp`, and `mpi-openmp` configure/build/test presets with separate build and install directories. The default build keeps MPI and OpenMP disabled.
- Export `ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` as interface targets with their own header sets and dependencies. Separate container storage, numerical algorithms and HDF5 adapters; foundations no longer inherit aggregate `ALPS::headers`. Numerical algorithms and public include names are preserved.
- Organize `src/alps/` by responsibility, with module-local headers, sources and tests. Configuration templates live in `cmake/config/`, and generated headers use `<build-dir>/generated/include/alps/`. Explicit header file sets preserve public include names. See the [module layout](src/alps/README.md).
- Export independently linkable runtime components `ALPS::utilities`, `ALPS::hdf5`, `ALPS::params`, `ALPS::osiris`, `ALPS::xml` and `ALPS::cli`. `ALPS::alps` links them transitively; Python packages carry one copy of each component. Params text/XML and older `Parameters` conversion adapters remain in `ALPS::alps`, with unchanged public header names. Archive formats and existing parsing behavior are preserved. Rebuild downstream binaries after the library splits.
- MaxEnt's solver and executable link the foundation components without `ALPS::alps`. The executable uses `ALPS::cli` and reads typed params directly from HDF5. Scientific calculations and stop-callback behavior are preserved.
- Group `src/tools/` commands by responsibility. Keep historical inactive sources without enabling them; executable names and installation components are preserved. `pconfig` now links only utilities.
- Require CMake 3.27 or newer and an externally installed Boost 1.76 or newer with CMake packages. The SDK requires C++17/C11 compilers, HDF5's C library, and LP64 BLAS/LAPACK; bundled Boost builds and alternate numerical integer/symbol ABIs are no longer supported.
- Export CMake targets for the installed SDK, applications and solver libraries. Downstream C++ projects link `ALPS::alps`; Python extensions sharing pyalps objects use `pyalps::runtime`. MaxEnt, CT-HYB and CT-INT Python wrappers link the SDK's solver libraries instead of compiling their implementations again.
- Make MPI opt-in with `ALPS_ENABLE_MPI=ON`. Standalone builds enable applications and native tests by default; embedded `add_subdirectory` builds default to the library alone. The default SDK uses shared libraries, as required by the Python bindings.
- Top-level single-configuration builds still default to `CMAKE_BUILD_TYPE=Release`; embedded builds leave the build type to the parent project.
- Support editable Python development against an installed SDK, with packaged native runtime resources and explicit downstream ABI checks. Redistributable wheels bundle native dependencies through platform repair and are tested after installation on fresh runners.
- Build the C++ and Fortran examples in `tutorials/examples/` as standalone projects against an installed SDK. Tutorial sources are installed only when the `tutorials` component is requested.
- Adopt ALPSCore's `<alps/utilities/...>` include spellings for `cast`, `stacktrace`, `signal`, `stringify`, `short_print`, `remove_cvr`, `type_wrapper`, `size`, `temporary_filename` and `make_copy`. `alps::ngs::signal` and `alps::ngs::stacktrace()` become `alps::signal` and `alps::stacktrace()`, and `ALPS_NGS_STRINGIFY` becomes `ALPS_STRINGIFY`. `alps::signal` moves from `ALPS::hdf5` to `ALPS::utilities`; open archives still close on a segmentation fault. The old spellings and names remain, deprecated, while `ALPS_LEGACY_API` is on (the default); see [the legacy layer](legacy/README.md). Rebuild downstream binaries.
- Retire the `alps/utility/` prefix. Its remaining shared helpers move to `<alps/utilities/...>`; `numeric_cast.hpp` moves to `<alps/numeric/numeric_cast.hpp>` and `cli.hpp` to `<alps/cli.hpp>`. `assign.hpp` and `set_zero.hpp` keep their spelling and now belong to the ALEA headers that use them. The unused `bitops.hpp` and `generator.hpp` are deprecated and installed only with the legacy layer. `factory.hpp` is removed: it has not compiled since before 3.0.0. See the [utilities status table](src/alps/utilities/README.md).
- Add ALPSCore's MPI wrapper `<alps/utilities/mpi.hpp>` (`alps::mpi::communicator`, `environment`, `broadcast`, `all_reduce`, `all_gather`) to `ALPS::utilities` in MPI builds, with its tests. It shares the namespace `alps::mpi` with `<alps/ngs/mpi.hpp>` without overlapping, and both can be used in one program.
- Reject HDF5 archive modes containing letters other than `r`, `w`, `a`, `c`, `l` and `m` with `alps::hdf5::wrong_mode`. Unknown letters, such as a mistyped `"wx"`, used to be ignored. The integer-property constructor likewise rejects bits other than `READ`, `WRITE`, `REPLACE`, `COMPRESS`, `LARGE` and `MEMORY`.
- Default-construct `alps::hdf5::archive` and attach a file later with `open(path, mode)`, also after `close()`; opening an open archive throws `alps::hdf5::archive_opened`. Copy assignment is deleted: it shared the file without counting the reference, so closing either archive could close the other's file. Copy construction is unchanged.
- Add `alps::hdf5::is_native_type<T>` and HDF5 storage for `boost::optional` in `<alps/hdf5/boost_optional.hpp>`, both from ALPSCore. An empty optional stores a default value with the attribute `alps_hdf5_optional_empty`, as ALPSCore does.
- Add `ALPS::aleav2`, ALEA v2: the ALEA library of the ALPS Core Libraries, copied from ALPSCore with its tests. It keeps ALPSCore's include spellings and namespace (`<alps/alea.hpp>`, `<alps/alea/*.hpp>`, `alps::alea`). It requires Eigen 3.3 or newer and is built only when Eigen is found; installed SDK consumers request it with `find_package(ALPS COMPONENTS aleav2)`. In MPI builds it includes ALPSCore's `alps::alea::mpi_reducer`. ALPS's own ALEA (`<alps/alea.h>`) remains the default: it is unchanged and not deprecated, and both can be used in the same program. Nothing in ALPS uses ALEA v2 yet; the two will be swapped only after a conformance suite passes on both. The citation catalog's new `aleav2` component credits the ALPS Core Libraries paper (Wallerberger et al., arXiv:1811.08331). See [its README](src/alps/aleav2/README.md) for provenance and license.
- Adopt ALPSCore's `<alps/numeric/vector_functions.hpp>` and add its `<alps/numeric/inf.hpp>`. Vector arithmetic in `alps::numeric` now treats a default-constructed vector as the zero vector, and `+=`, `-=`, `*=` and `/=` throw `std::runtime_error` on a size mismatch. Scaling by a scalar of another arithmetic type, unary `+` and `operator<<` keep working. The functors of `<alps/numeric/functional.hpp>` now come from `vector_functions.hpp`. Code that declares both `boost::numeric::operators::operator*` and `alps::numeric::operator*` (or `+`, `-`, `/`) in one scope must drop the Boost declaration. See [the numerics README](src/alps/numerics/README.md).
- Add `ALPS::accumulators`, the accumulators library of the ALPS Core Libraries, copied from ALPSCore with its tests (`alps::accumulators`, `<alps/accumulators.hpp>`, `<alps/accumulators/*.hpp>`). It is always built, is part of `ALPS::alps`, and installed SDK consumers can request it with `find_package(ALPS COMPONENTS accumulators)`. The NGS accumulators (`<alps/ngs/accumulator/*.hpp>`) are unchanged. `<alps/config.hpp>`, ALPSCore's name for the configuration header, now ships alongside `<alps/config.h>`, and ALPSCore's `<alps/numeric/{check_size,scalar,rectangularize,boost_array_functions}.hpp>` and, in MPI builds, `<alps/utilities/boost_mpi.hpp>` are added. The citation catalog's new `accumulators` component credits the ALPS Core Libraries paper. See [its README](src/alps/accumulators/README.md) for provenance and license.
- Add ALPSCore's parameters to `ALPS::params`, copied from ALPSCore with their tests: `alps::params_ns::params` and `alps::params_ns::dictionary` in `<alps/params.hpp>`, `<alps/dictionary.hpp>` and `<alps/params/*.hpp>`, with the iniparser library they read INI files with. Like the NGS `alps::params` of `<alps/ngs/params.hpp>`, which is unchanged, ALPSCore's header names its class `alps::params`, so one translation unit can include only one of the two headers. The citation catalog's new `params` component credits the ALPS Core Libraries paper. See [its README](src/alps/params/README.md) for provenance and licenses.
- `pyalps.loadBinningAnalysis` reads the binning analysis of results written by ALPSCore's accumulators (`error_bins`), as well as ALEA's `timeseries/logbinning` data. `pyalps.loadMeasurements` already read their means, errors and counts, and their parameters. A frozen file written by ALPSCore checks both. The returned `hdf5_path` now includes the `/` before the observable name.
- The hybridization solver (`hybridization`, `pyalps.cthyb`) gathers its parameter defaults in `hybridization::define_parameters`, so the `/parameters` it writes now list the defaults it used. A measurement that needs `N_MATSUBARA`, `N_LEGENDRE`, `N_nn`, `N_W` or `N_w2` now rejects the value 0 as well as a missing value, and an empty `U_MATRIX`, `MU_VECTOR` or `RET_INT_K` counts as unset.
- The interaction-expansion solver (`interaction`, `pyalps.ctint`) gathers its parameter defaults in `InteractionExpansion::define_parameters`, and its simulation and evaluation read the sizes with older spellings or computed defaults (`N_MATSUBARA` or `NMATSUBARA`, `N_TAU` or `N`, `FLAVORS` or `N_ORBITALS`, `NSELF`, `NMATSUBARA_MEASUREMENTS`, `MEASUREMENT_PERIOD`, `CONVERGENCE_CHECK_PERIOD`) through one set of functions. `ATOMIC` and `TAU_DISCRETIZATION_FOR_EXP` now act by value: `ATOMIC=0` no longer selects the atomic limit, and `TAU_DISCRETIZATION_FOR_EXP=0` computes the exponentials on the fly.
- MaxEnt gathers its parameter defaults in `MaxEntSimulation::define_parameters`. An empty `DATA` or `COVARIANCE_MATRIX` counts as unset.

### Migration

Public and vendored headers are retained. Tutorial consumers use the exported SDK targets. The ALPS 2.x consumer interfaces `UseALPS.cmake`, `UsePyALPS.cmake` and `include.mk` remain installed as deprecated compatibility layers over those targets. Configure downstream projects using [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code).

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
| `UseALPS.cmake`, `include.mk` | Imported SDK targets and explicit installation paths; add the installed `bin` directory to `PATH`. Both remain installed, deprecated: `${ALPS_LIBRARIES}` is `ALPS::alps`, and `UseALPS.cmake` no longer sets compilers, global flags or cache entries |
| `UsePyALPS.cmake` | [pyalps downstream CMake package](python/pyalps/README.md#downstream-native-extensions); the deprecated `alps_target_link_pyalps` links `pyalps::runtime` |
| `<alps/ngs/detail/export_sim_to_python.hpp>` | `<pyalps/export_simulation.hpp>`; the old header forwards to it |
| `<alps/ngs/{cast,stacktrace,signal,stringify,short_print}.hpp>`, `<alps/ngs/detail/{remove_cvr,type_wrapper}.hpp>`, `<alps/utility/{size,temporary_filename,make_copy}.hpp>` | `<alps/utilities/...>` with the same file name; the old headers forward to them while `ALPS_LEGACY_API` is on |
| `<alps/utility/{copyright,citations,os,vmusage,vectorio,encode,resize,data}.hpp>` | `<alps/utilities/...>` with the same file name; the old headers forward to them while `ALPS_LEGACY_API` is on |
| `<alps/utility/numeric_cast.hpp>`, `<alps/utility/cli.hpp>` | `<alps/numeric/numeric_cast.hpp>`, `<alps/cli.hpp>`; the old headers forward to them while `ALPS_LEGACY_API` is on |
| `<alps/utility/{bitops,generator}.hpp>` | No replacement; installed, deprecated, while `ALPS_LEGACY_API` is on. Use C++20 `<bit>` instead of `bitops` |
| `<alps/utility/factory.hpp>` | Removed. It did not compile in 3.0.0, so no code can depend on it |
| `<alps/ngs/parapack/*.h>` | Removed with their five sources: a parapack front end on the NGS parameters, which nothing used. The parapack scheduler `<alps/parapack/*.h>` is unchanged |
| `<alps/ngs/scheduler/proto/*.hpp>` | Removed: an unused prototype scheduler on the NGS parameters and `mcbase`. Use `alps::mcbase` with `<alps/mcmpiadapter.hpp>` |
| `alps::ngs::signal`, `alps::ngs::stacktrace()`, `ALPS_NGS_STRINGIFY` | `alps::signal`, `alps::stacktrace()`, `ALPS_STRINGIFY`; deprecated aliases remain while `ALPS_LEGACY_API` is on |
| Historical XML shell tools | Retained alongside `alps-xml`; see [XML tools](CONTRIBUTING.md#xml-resources-and-tools) |
| Top-level `import mpi` compatibility module | `import pyalps.mpi`; install the `mpi` extra for mpi4py |
| `applications/`, `tool/`, `test/`, `example/` | `src/apps/`, `src/tools/`, `tests/`, `tutorials/examples/` |
| Loose subsystem trees and transitional `src/alps/{common,runtime}/` | Semantic `src/alps/<module>/{include,src,tests}` ownership; see the [module map](src/alps/README.md#source-ownership) |
| `src/ietl/` | `src/alps/ietl/include/ietl/`; public `<ietl/...>` includes are unchanged |
| `src/alps/config.h.in`, `src/alps/version.h.in`, transitional `src/alps/common/config/` | `cmake/config/` configuration templates |
| Build-tree `src/alps/` generated headers | `<build-dir>/generated/include/alps/`; consumers use exported CMake targets |
| XML definitions and stylesheets in `lib/xml/` | `src/alps/resources/` in the source tree; installed under `share/alps/xml/` |
| `tutorials/ngs/5_export_python/` | `python/pyalps/examples/ising/` |
| Root `ALPS_VERSION.txt` | `cmake/ALPS_VERSION.txt`, still shared by the SDK and Python package |
