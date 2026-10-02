# Module layout and ALPSCore reconciliation

This layout prepares ALPS for incremental reconciliation with ALPSCore after the CMake modernization. It establishes ownership of existing code and tests. It does not import ALPSCore implementations, change scientific algorithms or claim that the modules can already be linked independently.

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

Each library module contributes sources to the existing `alps` target and a named `HEADERS` file set to `alps_headers`. The file sets supply build include roots and install the headers at their existing `<alps/...>` paths. Public template implementation headers under `include/alps/ngs/detail/` remain installed because consumers need them. Module `src/`, `tests/` and the physical nesting of `include/` are excluded from the SDK.

`ALPS::alps` remains the single runtime library, shared by default and required to be shared for Python extensions. `ALPS::headers` remains its compile interface, and `ALPS::maxent` remains the solver library with its public entry point declared in `<alps/solvers.hpp>`. Do not add component aliases that imply independently linkable libraries before those boundaries exist. Module CMake files are subdirectories of the main build, not standalone projects.

## Dependencies to reconcile

| Area | Current coupling | Work needed before independent libraries |
| --- | --- | --- |
| Utilities | Generated ALPS configuration/version headers; NGS configuration; `utility/vectorio.hpp` uses the parser and element traits | Define the small common configuration interface and separate parser-facing helpers from a foundational utility library |
| HDF5 | Utility casts/stack traces, shared configuration, numeric-container adapters; archive and NGS signal handler call each other | Keep `signal.cpp` and `ngs/signal.hpp` with HDF5 for now; separate generic signal handling from archive cleanup, then define adapter dependencies |
| NGS params | HDF5 serialization, utility helpers, numeric helpers; constructors and XML conversion use the older parameter/parser subsystem | Agree on typed access, missing/default values, iteration and input semantics; isolate old-format adapters before replacing an implementation |
| MaxEnt | NGS/mcbase, scheduler, parameters, observables, HDF5 and numerical libraries | Reduce framework coupling around the solver and compare input/output and numerical behavior against the candidate ALPSCore-based implementation |

Directory separation alone does not remove these dependencies. This change deliberately keeps the existing target dependency graph and one runtime/export configuration while making the next changes local and reviewable.

## Reconciliation sequence

1. **Establish this structural baseline.** Build the runtime and applications, run existing tests, install the SDK and build external consumers. Verify that public header paths, test identities and scientific implementations survive the moves.
2. **Compare HDF5 and params contracts.** Record namespace/symbol overlap, dependency and license requirements, archive formats, scalar/container conversions, parameter parsing and error behavior. Add compatibility fixtures before selecting implementations. Do not link overlapping ALPS/ALPSCore definitions into one process without an explicit symbol-ownership plan.
3. **Extract actual library boundaries.** Resolve the couplings above, then introduce real component targets with explicit dependencies and exports. Exercise each target from an installed consumer. Keep the Python runtime identity coherent across extensions.
4. **Use MaxEnt as the first scientific integration.** Route its parameter and archive access through the agreed interfaces. Retain the current numerical regression and add representative scientific/reference data comparisons with agreed tolerances before switching solvers.
5. **Extend the same pattern only as needed.** Reconcile the next subsystem after the initial contracts are stable; avoid moving all remaining directories merely for visual consistency.

Shared build policy stays in the top-level CMake files and `cmake/`. Module CMake files declare sources, public headers and local tests without a second project/version/dependency-discovery framework.

## Validation

`ALPS_BUILD_TESTING` controls both module-local and central tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its test. The existing test names are preserved; labels permit focused runs after building:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|maxent)$'
```

Run the full native suite for changes to these shared interfaces. The installed SDK consumer in `tests/cmake/consumer/` exercises HDF5, params and utilities through exported targets and rejects private module directories in the installation. The SDK contract suite also checks C++17/C++20 consumers, embedding and relocation. Python packaging continues to consume the installed SDK; validate its bindings before changing binary or API ownership.
