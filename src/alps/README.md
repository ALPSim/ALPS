# Module layout and ALPSCore reconciliation

This layout prepares ALPS for incremental reconciliation with ALPSCore after the CMake modernization. It establishes ownership of existing code and tests and extracts utilities, HDF5 and typed params as independently linkable libraries. It does not import ALPSCore implementations or change scientific algorithms.

## Source ownership

```text
src/alps/
  CMakeLists.txt          # ALPS::alps runtime, common compile interface and install rules
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
  parameter/             # Existing alps::Parameters implementation; not NGS params
  ...                    # Remaining ALPS subsystems retain their current layout
src/apps/maxent/
  CMakeLists.txt         # ALPS::maxent solver and maxent executable
  src/                   # Solver implementation and private headers
  cli/                   # Command-line entry point
  tests/                 # Numerical regression
tests/integration/hdf5/ # Archive compatibility with older parameters and observables
tests/integration/params/ # Text/XML adapter contracts across the runtime boundary
```

The module name `params` refers to `alps::params` (`<alps/ngs/params.hpp>`). The older `alps::Parameters` API (`<alps/parameter.h>`) remains in `parameter/`. Keep that distinction explicit when comparing ALPSCore APIs.

Utilities, HDF5 and typed params build the `alps_utilities`, `alps_hdf5` and `alps_params` libraries. The sources in `params/adapters/` contribute to `alps`, alongside the older parameter/parser implementation. Each module contributes a named `HEADERS` file set to the shared SDK compile interface `alps_headers`. The file sets supply build include roots and install the headers at their existing `<alps/...>` paths. Public template implementation headers under `include/alps/ngs/detail/` remain installed because consumers need them. Module `src/`, `tests/` and the physical nesting of `include/` are excluded from the SDK.

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
| MaxEnt | NGS/mcbase, scheduler, parameters, observables, HDF5 and numerical libraries | Reduce framework coupling around the solver and compare input/output and numerical behavior against the candidate ALPSCore-based implementation |

Directory separation alone does not remove these dependencies. Utilities, HDF5 and typed params now have tested binary boundaries. The params adapters keep the older parameter/parser dependency outside the typed runtime.

## Reconciliation sequence

1. **Establish this structural baseline.** Build the runtime and applications, run existing tests, install the SDK and build external consumers. Verify that public header paths, test identities and scientific implementations survive the moves.
2. **Compare HDF5 and params contracts.** Record namespace/symbol overlap, dependency and license requirements, archive formats, scalar/container conversions, parameter parsing and error behavior. Add compatibility fixtures before selecting implementations. Do not link overlapping ALPS/ALPSCore definitions into one process without an explicit symbol-ownership plan.
3. **Extend actual library boundaries.** Utilities, HDF5 and typed params are extracted and tested through their own targets. Keep adapter ownership explicit when reconciling implementations. Exercise each target from an installed consumer. Keep the Python runtime identity coherent across extensions.
4. **Use MaxEnt as the first scientific integration.** Route its parameter and archive access through the agreed interfaces. Retain the current numerical regression and add representative scientific/reference data comparisons with agreed tolerances before switching solvers.
5. **Extend the same pattern only as needed.** Reconcile the next subsystem after the initial contracts are stable; avoid moving all remaining directories merely for visual consistency.

Shared build policy stays in the top-level CMake files and `cmake/`. Module CMake files declare sources, public headers and local tests without a second project/version/dependency-discovery framework.

## Validation

`ALPS_BUILD_TESTING` controls both module-local and central tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its test. The existing test names are preserved; labels permit focused runs after building:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|maxent)$'
```

Run the full native suite for changes to these shared interfaces. The installed SDK consumer in `tests/cmake/consumer/` exercises HDF5 and params through `ALPS::alps`, and separate executables link only `ALPS::utilities`, `ALPS::hdf5` or `ALPS::params`. HDF5 and params module tests also link only their components; integration tests for older parameters, observables and input adapters keep `ALPS::alps`. The adapter contract also runs through the installed SDK. A static embedding test builds and runs only the three extracted components without building the simulation runtime. The installed consumer also rejects private module directories in the installation. The SDK contract suite also checks C++17/C++20 consumers, embedding and relocation. Python packaging continues to consume the installed SDK; validate its bindings before changing binary or API ownership.
