# Module layout and ALPSCore reconciliation

ALPS sources are grouped by responsibility to prepare MaxEnt, HDF5 and typed params for ALPSCore reconciliation. Utilities, HDF5, params, Osiris, XML and command-line parsing have separate runtime libraries. Other source modules contribute to the existing `ALPS::alps` library or its shared compile interface. This cleanup imports no ALPSCore implementation and preserves scientific algorithms.

The [reconciliation report](../../doc/ALPSCore-reconciliation.md) pins the ALPSCore reference, records measured archive/params compatibility and defines migration gates. Reproduce that comparison with the [separate-process probes](../../tests/reconciliation/README.md).

## Source ownership

Each module uses `include/`, `src/` and `tests/` where applicable. Public include spellings describe the API, independently of the physical owner: for example, NGS measurement headers live in `alea/include/alps/ngs/`, while typed parameters live in `params/include/alps/ngs/`.

| Source module under `src/alps/` | Responsibility | Binary or compile owner |
| --- | --- | --- |
| `utilities/` | Utility functions, general helpers, type traits and NGS configuration helpers | `ALPS::utilities`, `ALPS::headers` |
| `containers/` | Fixed-capacity and ALPS multi-array containers | `ALPS::headers` |
| `numerics/` | Numerical helpers and matrix/vector interfaces | `ALPS::headers` |
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
| `solvers/` | Shared `<alps/solvers.hpp>` declarations for MaxEnt and CT-QMC | `ALPS::headers` |
| `resources/` | XML definitions and stylesheets | Installed data component |

The broad `common/` and `runtime/` source groups are removed. ALPS and IETL configuration templates live in `cmake/config/`; generated public headers live under `<build-dir>/generated/include/alps/`. MaxEnt remains in `src/apps/maxent/{src,cli,tests}`. Subsystem tests follow their source owner; cross-module integration, SDK, Python, CLI, packaging and reconciliation checks stay under root `tests/`. Relocation preserves registered test names and leaves inactive fixtures inactive.

All first-party public headers have explicit CMake `HEADERS` file sets. These declare build include roots and preserve installed `<alps/...>` and `<ietl/...>` paths; consumers use exported targets rather than broad source-tree include roots. Public template definitions remain installed, while private source files, tests and physical `include/` nesting do not enter the SDK.

## Library boundaries

`ALPS::utilities` owns utility symbols and links Boost.Filesystem and platform threads without the simulation runtime, HDF5 or BLAS/LAPACK.

`ALPS::hdf5` owns archive symbols, exception exports, shared archive state and the NGS signal handler that closes archives. It links utilities, HDF5, Boost.Filesystem, Boost.Thread and platform threads. Archive and signal code remain together to preserve cleanup behavior. A parallel HDF5 provider can bring its own MPI dependency.

`ALPS::params` owns typed values, proxies, checkpoint I/O and the `paramvalue_source` interface used by Python bindings. It links HDF5 and Boost.Serialization; MPI builds additionally link MPI and Boost.MPI. Its parameter-file constructor, XML input and older `Parameters` conversion remain in `params/adapters/`, compiled into `ALPS::alps`. The older `alps::Parameters` implementation itself belongs to `legacy_parameters/`.

`ALPS::osiris` owns dump/process APIs, XDR symbols and the communication state used by `comm_init()` and `is_master()`. It links Boost.Serialization/Filesystem and, when enabled, MPI. This gives communication state one owner and permits MaxEnt to preserve existing diagnostic gating without linking the simulation runtime.

`ALPS::xml` owns XML parsing, attributes, handlers, output streams and stylesheet lookup, with Boost.Filesystem and Boost.Regex dependencies. XML file-to-parameter conversion still belongs to the params adapters in `ALPS::alps`. The separate `plotting/` header owner keeps `<alps/plot.h>` and its older `Parameters` dependency outside the XML component.

`ALPS::cli` owns the existing `mcoptions` and `parseargs` implementations, linking utilities and Boost.ProgramOptions. Their installed headers remain `<alps/ngs/mcoptions.hpp>` and `<alps/parseargs.hpp>`. The two existing option grammars, defaults, filename rules and error behavior are preserved; these parsers do not read parameter files.

`ALPS::maxent` links params, HDF5, utilities, Osiris and its Boost/numerical providers. It owns its deterministic run loop; its numerical calculations and stop-callback ordering are preserved. The `maxent` executable adds `ALPS::cli` for the existing `mcoptions` parser and reads typed params from HDF5 directly. Neither the solver nor executable links `ALPS::alps`. Its public callable API remains `<alps/solvers.hpp>`.

`ALPS::alps` links the extracted runtime components publicly. These libraries follow `BUILD_SHARED_LIBS`, with component-specific generated export headers. Python extensions require shared runtime libraries and package one copy of every component in `ALPS_RUNTIME_TARGETS`. Rebuild downstream binaries after the XML and CLI extractions, as after the earlier library splits.

Physical ownership does not imply independent linkability. `ALPS::headers` still exposes the aggregate compile interface; numerical headers depend on parser support, and other include cycles remain. The semantic modules make these dependencies visible without claiming that each can already be configured or linked alone. Package discovery still checks the full SDK dependency set. Shared build policy remains in the root CMake files and `cmake/`.

## Architecture checks

CMake generates `alps-module-manifest.json` from module declarations, public-header file sets and source lists. The [architecture checker](../../.github/scripts/check_module_architecture.py) checks production-file ownership, public include ownership, declared include dependencies and public/private boundaries. It reports observed cycles, nonliteral includes and documented unresolved includes for review. It does not replace compilation, link-dependency checks or scientific validation.

After configuring the build, run:

```sh
python .github/scripts/check_module_architecture.py \
  --manifest _build/default/alps-module-manifest.json \
  --write-report _build/default/alps-module-report.json
```

Update the owning module's CMake declarations when adding files or dependencies; the manifest is generated rather than hand-maintained. A declared existing cycle is visible architectural debt, not evidence of independent components.

### Recorded architectural debt

After the XML and CLI extractions, `alps-module-architecture.json` inventories 26 source owners, 517 public include spellings and 626 production files. The owners include `cli`, `plotting` and separate MaxEnt solver/executable owners for `src/apps/maxent/src/` and `src/apps/maxent/cli/`. These are ownership counts, not counts of independent libraries or passing tests. The earlier code checkpoint `627d4f500` had 24 owners, before the CLI and plotting modules were separated.

The current observed include graph contains a five-module cycle (`containers`, `hdf5`, `numerics`, `utilities`, `xml`) and a two-module cycle (`expression`, `legacy_parameters`). The earlier eight-module cycle at `627d4f500` has therefore narrowed, but the remaining cycles still need deliberate reconciliation. Dependency declarations constrain new include edges. Regenerate the report after changing module ownership or dependencies.

The report inventories 82 exact unresolved file/include pairs already present at baseline `f27ed2316`: 80 in dormant accumulator code and two in the optional `USE_LATTICE_CONSTANT_2D` graph backend. Each exemption names its file, include and reason; they do not establish support for those inactive paths. Resolve or remove these dependencies deliberately rather than adding broad exclusions.

With `ALPS_BUILD_TESTING=ON`, CTest runs `module_architecture` and writes `<build-dir>/alps-module-architecture.json`; this requires a Python interpreter ≥ 3.10. Builds with testing disabled do not need Python for module configuration or manifest generation.

## Reconciliation sequence

1. **Validate the boundaries.** Build applications and native tests, install the SDK and exercise individual component consumers, including Osiris, XML, CLI and MaxEnt. Verify public headers, test identities, relocation and Python runtime ownership. Check that the MaxEnt executable's runtime dependencies exclude `ALPS::alps`.
2. **Resolve HDF5 and params contracts.** Use the pinned cross-read fixtures to address archive metadata, unsupported value decoding, conversions, exception behavior and resource lifetime. Keep older scientific input adapters explicit.
3. **Use MaxEnt as the first scientific workload.** Extend the existing linear-grid regression with representative kernels, grids, covariance and reference spectra before selecting replacement implementations. ALPSCore's pinned repository contains no MaxEnt solver to import.
4. **Replace one owner at a time.** Never load overlapping ALPS and ALPSCore archive implementations into one process. Preserve the agreed API, serialization and Python contracts, then rebuild dependents.
5. **Extract further libraries when dependencies permit.** Use the declared and observed module graph to choose the next boundary; source directories alone do not establish it.

## Validation

`ALPS_BUILD_TESTING` controls module-local and central native tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its regression. After building, focused existing labels include:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|osiris|parser|cli|maxent)$'
```

Run the full native suite for shared-interface changes. The SDK consumers in `tests/cmake/consumer/` exercise aggregate and component links; integration tests for older inputs and observables keep `ALPS::alps`. Check shared/static consumers, installed-header ownership, SDK relocation and Python extension interoperability. The CLI contract checks existing defaults, option spellings, filename rules, help and execution-mode handling. The MaxEnt regression also exercises early stop, callback exceptions and completion ordering. These are validation requirements; this ownership map does not assert that every platform/configuration has passed them.
