# Module layout and ALPSCore reconciliation

This layout prepares ALPS for incremental reconciliation with ALPSCore after the CMake modernization. Utilities, HDF5 and typed params have independently linkable libraries. The remaining runtime, shared headers and Fortran bridge also have explicit source, include and test locations. This filesystem organization preserves existing binary ownership and scientific implementations.

The [reconciliation report](../../doc/ALPSCore-reconciliation.md) pins an ALPSCore
revision, records measured archive/params compatibility, and defines migration
gates. Reproduce the comparison with the [separate-process probes](../../tests/reconciliation/README.md).

## Source ownership

```text
src/alps/
  CMakeLists.txt          # Shared compile interface, component coordination and exports
  common/
    CMakeLists.txt        # Shared public headers and generated-header inputs
    config/              # ALPS config/version and IETL configuration templates
    include/alps/        # Numeric/container traits, shared helpers, solvers.hpp
    include/ietl/        # IETL headers
    tests/               # numeric/ and fixed_capacity/
  utilities/
    CMakeLists.txt        # Implementation sources and explicit public-header file set
    include/alps/        # Existing utility/ and selected low-level ngs/ include paths
    src/
    tests/
  hdf5/
    CMakeLists.txt
    include/alps/        # hdf5.hpp, hdf5/ and ngs/signal.hpp
    src/
    tests/               # Includes reference archives and serialization test helpers
  params/
    CMakeLists.txt
    include/alps/ngs/    # NGS params and its implementation headers
    src/                 # Typed params, proxy and value implementation
    adapters/            # Text/XML and older Parameters bridges, owned by ALPS::alps
    tests/
  runtime/
    CMakeLists.txt        # Remaining ALPS::alps implementation
    include/alps/        # Existing subsystem headers and public template definitions
    src/                 # alea/, expression/, lattice/, model/, ngs/, osiris/, ...
    tests/               # Subsystem-native tests, preserving their existing names
  fortran/
    CMakeLists.txt        # Existing ALPS::fortran bridge
    include/alps/fortran/
    src/
  resources/             # XML definitions and stylesheets
src/apps/maxent/
  CMakeLists.txt         # ALPS::maxent solver and maxent executable
  src/                   # Solver implementation and private headers
  cli/                   # Command-line entry point
  tests/                 # Numerical regression
tests/integration/hdf5/ # Archive compatibility with older parameters and observables
tests/integration/params/ # Text/XML adapter contracts across the runtime boundary
tests/                   # SDK, Python, CLI, build-policy and reconciliation checks
```

The module name `params` refers to `alps::params` (`<alps/ngs/params.hpp>`). The older `alps::Parameters` API (`<alps/parameter.h>`) belongs to `runtime/include/alps/parameter/` and `runtime/src/parameter/`. Keep that distinction explicit when comparing ALPSCore APIs.

Utilities, HDF5 and typed params build the `alps_utilities`, `alps_hdf5` and `alps_params` libraries. The sources in `params/adapters/` contribute to `alps`, alongside the older parameter/parser implementation in `runtime/`. The Fortran bridge retains its own `ALPS::fortran` target. These library boundaries are unchanged by the filesystem pass.

Every first-party public header is listed explicitly in a CMake `HEADERS` file set on the shared SDK compile interface `alps_headers`. The file sets supply component include roots and preserve installed `<alps/...>` and `<ietl/...>` paths. Public template definitions remain installed even when they are implementation details. Generated headers live under `<build-dir>/generated/include/alps/`; their inputs live in `common/config/`. Builds use the declared include roots rather than broad source-tree or build-tree `src/` include paths. Module sources, tests, configuration templates and physical `include/` nesting do not enter the SDK.

`common/` owns shared source headers, not an independently linkable API or library. It holds numerical/container helpers, type traits, fixed-capacity containers, selected NGS configuration helpers and the solver declarations shared by MaxEnt and CT-QMC. Some numerical headers still include runtime parser headers. `ALPS::headers` exposes these compile dependencies together; moving headers does not remove their dependencies or require a simulation-runtime link.

`runtime/` preserves the existing subsystem subdivisions under `include/alps/`, `src/` and `tests/`. Its test areas are alea, graph, lattice, model, NGS, Osiris, older parameters, parapack, parser, random and accumulator. Common numerical and fixed-capacity tests live in `common/tests/`. Relocation preserves test registration and names; inactive test fixtures remain inactive. Cross-module integration and build/package/Python checks stay under root `tests/`.

`ALPS::utilities` owns utility symbols and has no link dependency on `ALPS::alps`, HDF5 or the numerical libraries. It links Boost.Filesystem and the platform threading library.

`ALPS::hdf5` owns archive symbols, exception exports, shared archive state and the NGS signal handler that closes archives. It links utilities, HDF5, Boost.Filesystem, Boost.Thread and platform threads, with no direct dependency on `ALPS::alps`, MPI or BLAS/LAPACK; a parallel HDF5 provider can bring its own MPI dependency. Archive and signal code remain together to preserve cleanup behavior; their mutual calls are internal to this component.

`ALPS::params` owns typed values, lookup, iteration, parameter proxies, HDF5 checkpoint I/O and the `paramvalue_source` interface used by bindings. It links HDF5 and Boost.Serialization; MPI-enabled builds additionally link MPI and Boost.MPI for broadcast. It has no link dependency on `ALPS::alps` or BLAS/LAPACK.

The parameter-file constructor `params(boost::filesystem::path const&)`, `make_parameters_from_xml` and `make_deprecated_parameters` use the older parameter/parser implementation and remain in `ALPS::alps`. Link that target when using these adapters. `params` methods have individual export annotations: typed operations use `ALPS_PARAMS_DECL`, while the file constructor uses `ALPS_DECL`. This permits the constructor to live in a different Windows DLL without a dependency cycle. Public source APIs and input semantics are preserved.

`ALPS::alps` links all three components publicly, so existing source consumers keep the same target. All four runtime libraries follow `BUILD_SHARED_LIBS`; Python extensions require shared libraries and package one copy of each component listed in `ALPS_RUNTIME_TARGETS`. Utilities, HDF5 and params use `ALPS_UTILITIES_DECL`, `ALPS_HDF5_DECL` and `ALPS_PARAMS_DECL` with generated export headers, including for Windows DLLs. Downstream binaries must be rebuilt after the splits.

`ALPS::headers` remains a shared compile interface, and `ALPS::maxent` remains the solver library with its public entry point declared in `<alps/solvers.hpp>`. Package discovery still checks the full SDK dependency set. Module CMake files are subdirectories of the main build, not standalone projects.

## Dependencies to reconcile

| Area | Current coupling | Next work |
| --- | --- | --- |
| Utilities | Independently linkable; shared SDK configuration and header-only numeric/container traits remain compile dependencies. The unused parser include in `vectorio.hpp` is removed | Narrow the shared compile interface and package dependency discovery if a separately configurable utility package is needed |
| HDF5 | Independently linkable; utility casts/stack traces, shared configuration and numeric-container adapters. Archive and NGS signal cleanup have one runtime owner | Compare archive contracts with ALPSCore; assess signal ownership and adapter dependencies before implementation replacement |
| NGS params | Independently linkable typed runtime; HDF5 and header-only numeric helpers. Text/XML and older `Parameters` conversion are isolated adapters in `ALPS::alps` | Compare typed access, missing/default values, iteration, value conversions and input semantics against ALPSCore; keep adapter contracts explicit before replacing an implementation |
| MaxEnt | Typed params already in use; mcbase execution, Osiris diagnostic gating, CLI mcoptions, HDF5 and BLAS/LAPACK | Separate execution/diagnostics/CLI dependencies and use the solver as a scientific acceptance workload for reconciled foundations; the pinned ALPSCore repo contains no MaxEnt implementation |

Directory separation alone does not remove these dependencies. Utilities, HDF5 and typed params have tested binary boundaries. The params adapters keep the older parameter/parser dependency outside the typed runtime. `common/` and the runtime subsystem folders describe source ownership; they do not introduce further library boundaries.

## Reconciliation sequence

1. **Establish this structural baseline.** Build the runtime and applications, run existing tests, install the SDK and build external consumers. Verify that public header paths, test identities and scientific implementations survive the moves.
2. **Compare HDF5 and params contracts.** Record namespace/symbol overlap, dependency and license requirements, archive formats, scalar/container conversions, parameter parsing and error behavior. Add compatibility fixtures before selecting implementations. Do not link overlapping ALPS/ALPSCore definitions into one process without an explicit symbol-ownership plan.
3. **Extend actual library boundaries.** Utilities, HDF5 and typed params are extracted and tested through their own targets. Keep adapter ownership explicit when reconciling implementations. Exercise each target from an installed consumer. Keep the Python runtime identity coherent across extensions.
4. **Use MaxEnt as the first scientific integration.** Route its parameter and archive access through the agreed interfaces. Retain the current numerical regression and add representative scientific/reference data comparisons with agreed tolerances before switching solvers.
5. **Extract further libraries only as needed.** Reconcile the next subsystem after the initial contracts are stable. The consistent filesystem layout makes its sources easier to locate without implying that each historical subsystem is independently linkable.

Shared build policy stays in the top-level CMake files and `cmake/`. Module CMake files declare sources, public headers and local tests without a second project/version/dependency-discovery framework.

## Validation

`ALPS_BUILD_TESTING` controls both module-local and central tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its test. The existing test names are preserved; labels permit focused runs after building:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|maxent)$'
```

Run the full native suite for changes to these shared interfaces. The installed SDK consumer in `tests/cmake/consumer/` exercises HDF5 and params through `ALPS::alps`, and separate executables link only `ALPS::utilities`, `ALPS::hdf5` or `ALPS::params`. HDF5 and params module tests also link only their components; integration tests for older parameters, observables and input adapters keep `ALPS::alps`. The adapter contract also runs through the installed SDK. A static embedding test builds and runs only the three extracted components without building the simulation runtime. The installed consumer also rejects private module directories in the installation. The SDK contract suite also checks C++17/C++20 consumers, embedding and relocation. Python packaging continues to consume the installed SDK; validate its bindings before changing binary or API ownership.
