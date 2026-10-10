# C++ module layout

Sources use `include/`, `src/` and `tests/` where applicable. Public include spellings describe the API independently of the physical owner: NGS measurement headers live in `alea/include/alps/ngs/`, while typed parameters live in `params/include/alps/ngs/`.

## Source ownership

| Source module under `src/alps/` | Responsibility | Binary or compile owner |
| --- | --- | --- |
| `utilities/` | Utility functions, general helpers, type traits and NGS configuration helpers | `ALPS::utilities` |
| `containers/` | Fixed-capacity containers, ALPS multi-array storage and Boost serialization | `ALPS::containers` |
| `numerics/` | Numerical helpers, array mathematics and matrix/vector interfaces | `ALPS::numerics` |
| `numeric_io/` | HDF5 adapters for numerical matrices and vectors | `ALPS::numeric_io` |
| `ietl/` | Iterative eigensolver headers under `include/ietl/` | `ALPS::headers` |
| `hdf5/` | Archive API, container adapters, shared context registry and signal cleanup | `ALPS::hdf5` |
| `params/` | Typed values, lookup/proxies, iteration and HDF5 checkpoints | `ALPS::params`; `adapters/` contributes to `ALPS::alps` |
| `osiris/` | Dump serialization, process/communication state and XDR implementation | `ALPS::osiris` |
| `xml/` | XML parsing, handlers, attributes and output streams | `ALPS::xml` |
| `cli/` | Existing `mcoptions` and `parseargs` command-line grammars | `ALPS::cli` |
| `plotting/` | `<alps/plot.h>` output helpers combining XML and older parameters | `ALPS::headers` |
| `legacy_parameters/`, `expression/` | Older `alps::Parameters` and expression evaluation | `ALPS::alps` |
| `graph/`, `lattice/`, `model/` | Graph helpers, lattice definitions and physical models | `ALPS::headers`, `ALPS::alps` |
| `random/` | Random generators and their factories | `ALPS::alps` |
| `alea/`, `accumulators/` | Observable/result facilities and accumulator implementations | `ALPS::alps` |
| `mc/`, `scheduler/`, `parapack/` | Simulation API, execution and scheduling | `ALPS::alps` |
| `fortran/` | C++ bridge with public headers in `include/alps/fortran/` | `ALPS::fortran` |
| `solvers/` | Shared `<alps/solvers.hpp>` declarations for MaxEnt and CT-QMC | `ALPS::solver_headers` |
| `resources/` | XML definitions and stylesheets | Installed data component |

Configuration templates live in `cmake/config/`; generated headers live under `<build-dir>/generated/include/alps/`. Each owner's explicit CMake `HEADERS` file set preserves installed `<alps/...>` and `<ietl/...>` include paths. Public template definitions are installed; private sources and tests are not.

Subsystem tests follow their source owner. Cross-module integration, SDK, Python, CLI and packaging tests stay under root `tests/`. MaxEnt's implementation, CLI and tests live in `src/apps/maxent/{src,cli,tests}`. Tool groups are described in [src/tools/README.md](../tools/README.md).

## Library boundaries

`ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` are interface targets. `ALPS::containers` supplies storage, while `ALPS::numerics` adds numerical algorithms and BLAS/LAPACK. Numerical archive adapters belong to `ALPS::numeric_io`, which combines numerics and HDF5. See the [SDK usage and migration instructions](../../CONTRIBUTING.md#consuming-the-c-sdk).

Utilities, HDF5, typed params, Osiris, XML and command-line parsing are separate runtime libraries. `ALPS::alps` links them transitively. The params text/XML and older `Parameters` conversion adapters remain in `params/adapters/` and are compiled into `ALPS::alps`; their unchanged public header names are exposed through the aggregate interface.

`ALPS::maxent` links the foundations, Osiris and numerical providers. Its executable adds `ALPS::cli` and reads HDF5 params directly. Both use the shared `<alps/solvers.hpp>` declarations and link without `ALPS::alps`.

The remaining simulation modules contribute to `ALPS::alps` or its aggregate compile interface, `ALPS::headers`; physical ownership does not make each directory an independent library. Foundation targets declare their own dependencies without inheriting that aggregate interface. Package discovery currently checks the complete SDK dependency set.

Runtime libraries follow `BUILD_SHARED_LIBS` and have component-specific symbol exports. Python bindings require shared libraries and package one copy of each runtime component. Rebuild downstream binaries after changing the SDK or its dependency stack.

## Validation

Run the native suite for shared-interface changes. Component labels support focused checks:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|osiris|parser|cli|maxent)$'
```

The installed-SDK consumers in `tests/cmake/` check aggregate/component links, shared/static builds, header ownership, relocation and downstream extension interoperability. CLI and MaxEnt regressions check existing input behavior and scientific results. See [CONTRIBUTING.md](../../CONTRIBUTING.md#run-the-tests) for the complete development workflow.

Legacy text parameter loading is explicit: include `<alps/ngs/params_from_file.hpp>` and call `alps::params_from_file(path)` while linking `ALPS::alps`. The foundation `ALPS::params` supplies typed values and HDF5 persistence without the legacy text parser.
