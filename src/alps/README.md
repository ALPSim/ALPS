# Module layout and ALPSCore reconciliation

This layout prepares ALPS for incremental reconciliation with ALPSCore after the CMake modernization. It establishes ownership of existing code and tests and extracts utilities and HDF5 as independently linkable libraries. It does not import ALPSCore implementations or change scientific algorithms.

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
    src/
    tests/
  parameter/             # Existing alps::Parameters implementation; not NGS params
  ...                    # Remaining ALPS subsystems retain their current layout
src/apps/maxent/
  CMakeLists.txt         # ALPS::maxent solver and maxent executable
  src/                   # Solver implementation and private headers
  cli/                   # Command-line entry point
  tests/                 # Numerical regression
tests/integration/hdf5/ # Archive compatibility with older parameters and observables
```

The module name `params` refers to `alps::params` (`<alps/ngs/params.hpp>`). The older `alps::Parameters` API (`<alps/parameter.h>`) remains in `parameter/`. Keep that distinction explicit when comparing ALPSCore APIs.

Utilities and HDF5 build the `alps_utilities` and `alps_hdf5` libraries; params still contributes sources to `alps`. Each module contributes a named `HEADERS` file set to the shared SDK compile interface `alps_headers`. The file sets supply build include roots and install the headers at their existing `<alps/...>` paths. Public template implementation headers under `include/alps/ngs/detail/` remain installed because consumers need them. Module `src/`, `tests/` and the physical nesting of `include/` are excluded from the SDK.

`ALPS::utilities` owns utility symbols and has no link dependency on `ALPS::alps`, HDF5 or the numerical libraries. It links Boost.Filesystem and the platform threading library.

`ALPS::hdf5` owns archive symbols, exception exports, shared archive state and the NGS signal handler that closes archives. It links utilities, HDF5, Boost.Filesystem, Boost.Thread and platform threads, with no direct dependency on `ALPS::alps`, MPI or BLAS/LAPACK; a parallel HDF5 provider can bring its own MPI dependency. Archive and signal code remain together to preserve cleanup behavior; their mutual calls are internal to this component.

`ALPS::alps` links both components publicly, so existing source consumers keep the same target. All three follow `BUILD_SHARED_LIBS`; Python extensions require shared libraries and package one copy of each component listed in `ALPS_RUNTIME_TARGETS`. Utilities and HDF5 use `ALPS_UTILITIES_DECL` and `ALPS_HDF5_DECL` with generated export headers, including for Windows DLLs. Downstream binaries must be rebuilt after the splits.

`ALPS::headers` remains a shared compile interface, and `ALPS::maxent` remains the solver library with its public entry point declared in `<alps/solvers.hpp>`. Package discovery still checks the full SDK dependency set. Module CMake files are subdirectories of the main build, not standalone projects; params does not yet have an independent library target.

## Dependencies to reconcile

| Area | Current coupling | Next work |
| --- | --- | --- |
| Utilities | Independently linkable; shared SDK configuration and header-only numeric/container traits remain compile dependencies. The unused parser include in `vectorio.hpp` is removed | Narrow the shared compile interface and package dependency discovery if a separately configurable utility package is needed |
| HDF5 | Independently linkable; utility casts/stack traces, shared configuration and numeric-container adapters. Archive and NGS signal cleanup have one runtime owner | Compare archive contracts with ALPSCore; assess signal ownership and adapter dependencies before implementation replacement |
| NGS params | HDF5 serialization, utility helpers, numeric helpers; constructors and XML conversion use the older parameter/parser subsystem | Agree on typed access, missing/default values, iteration and input semantics; isolate old-format adapters before replacing an implementation |
| MaxEnt | NGS/mcbase, scheduler, parameters, observables, HDF5 and numerical libraries | Reduce framework coupling around the solver and compare input/output and numerical behavior against the candidate ALPSCore-based implementation |

Directory separation alone does not remove these dependencies. Utilities and HDF5 now have tested binary boundaries. Params retains its existing runtime ownership until its older parameter/parser coupling is resolved.

## Reconciliation sequence

1. **Establish this structural baseline.** Build the runtime and applications, run existing tests, install the SDK and build external consumers. Verify that public header paths, test identities and scientific implementations survive the moves.
2. **Compare HDF5 and params contracts.** Record namespace/symbol overlap, dependency and license requirements, archive formats, scalar/container conversions, parameter parsing and error behavior. Add compatibility fixtures before selecting implementations. Do not link overlapping ALPS/ALPSCore definitions into one process without an explicit symbol-ownership plan.
3. **Extend actual library boundaries.** Utilities and HDF5 are extracted and tested through their own targets. Resolve the params coupling above, then introduce its component target with explicit dependencies and exports. Exercise each target from an installed consumer. Keep the Python runtime identity coherent across extensions.
4. **Use MaxEnt as the first scientific integration.** Route its parameter and archive access through the agreed interfaces. Retain the current numerical regression and add representative scientific/reference data comparisons with agreed tolerances before switching solvers.
5. **Extend the same pattern only as needed.** Reconcile the next subsystem after the initial contracts are stable; avoid moving all remaining directories merely for visual consistency.

Shared build policy stays in the top-level CMake files and `cmake/`. Module CMake files declare sources, public headers and local tests without a second project/version/dependency-discovery framework.

## Validation

`ALPS_BUILD_TESTING` controls both module-local and central tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its test. The existing test names are preserved; labels permit focused runs after building:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|maxent)$'
```

Run the full native suite for changes to these shared interfaces. The installed SDK consumer in `tests/cmake/consumer/` exercises HDF5 and params through `ALPS::alps`, and separate executables link only `ALPS::utilities` or `ALPS::hdf5`. HDF5 module tests also link only its component; integration tests for parameters and observables keep `ALPS::alps`. A static embedding test builds and runs only the two extracted components without building the simulation runtime. The installed consumer also rejects private module directories in the installation. The SDK contract suite also checks C++17/C++20 consumers, embedding and relocation. Python packaging continues to consume the installed SDK; validate its bindings before changing binary or API ownership.
