# Module layout and ALPSCore reconciliation

ALPS sources are grouped by responsibility to prepare MaxEnt, HDF5 and typed params for ALPSCore reconciliation. Utilities, HDF5, params, Osiris, XML and command-line parsing have separate runtime libraries. Other source modules contribute to the existing `ALPS::alps` library or its shared compile interface. This cleanup imports no ALPSCore implementation and preserves scientific algorithms.

The [separate-process probes](../../tests/reconciliation/README.md) record the pinned ALPSCore reference and reproduce measured archive/params compatibility checks.

Before replacing HDF5, agree on the generic serialization contract with the new ALEA maintainer: buffer ownership and lifetime, shapes/types, error behavior and the owner of the HDF5 implementation. Start one measurement-application pilot in parallel with params reconciliation, and explicitly list the legacy result and checkpoint formats it must preserve. MaxEnt remains an acceptance test for the foundations; it does not exercise ALEA migration. `src/alps/alea/` still contains the legacy `Observable`/`ObservableSet` implementation; this pass imports neither new ALEA nor a compatibility shim. The pinned ALPSCore probes disable ALEA, so they do not validate ALEA compatibility.

## Source ownership

Each module uses `include/`, `src/` and `tests/` where applicable. Public include spellings describe the API, independently of the physical owner: for example, NGS measurement headers live in `alea/include/alps/ngs/`, while typed parameters live in `params/include/alps/ngs/`.

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

The broad `common/` and `runtime/` source groups are removed. ALPS and IETL configuration templates live in `cmake/config/`; generated public headers live under `<build-dir>/generated/include/alps/`. `ALPS::configuration` carries the shared configuration headers and compile requirements. MaxEnt remains in `src/apps/maxent/{src,cli,tests}`. Subsystem tests follow their source owner; cross-module integration, SDK, Python, CLI, packaging and reconciliation checks stay under root `tests/`. Relocation preserves registered test names and leaves inactive fixtures inactive.

All first-party public headers have explicit CMake `HEADERS` file sets. These declare build include roots and preserve installed `<alps/...>` and `<ietl/...>` paths; consumers use exported targets rather than broad source-tree include roots. Public template definitions remain installed, while private source files, tests and physical `include/` nesting do not enter the SDK.

## Library boundaries

`ALPS::configuration`, `ALPS::containers`, `ALPS::numerics`, `ALPS::numeric_io` and `ALPS::solver_headers` are exported interface targets, not new binary libraries. Foundation targets own their public header sets and declare their actual dependencies; they do not inherit `ALPS::headers`. The aggregate compile interface instead collects these foundations and the remaining simulation headers.

`ALPS::containers` provides storage without the mathematical umbrella. Array functions/operators, `<alps/multi_array.hpp>`, `<alps/functional.h>` and `<alps/utility/numeric_cast.hpp>` now belong to numerics, with their implementations and public include names preserved. Utility array resizing, abbreviated printing and MPI helpers include only the container facilities they use. `ALPS::numerics` adds utilities, containers and BLAS/LAPACK without HDF5 or XML.

`ALPS::numeric_io` combines numerics and HDF5. Include `<alps/hdf5/matrix.hpp>` or `<alps/hdf5/numeric_vector.hpp>` explicitly and link `ALPS::numeric_io` when using these archive adapters; `<alps/numeric/matrix.hpp>` no longer includes the matrix archive adapter. Diagonal/deprecated matrix/vector HDF5 `save`/`load` members and matrix `write_xml`/XML insertion overloads, unused in this repository, are removed without replacement adapters. This pass therefore includes API removals as well as ownership moves; it changes no numerical algorithms and imports no ALPSCore implementation.

`ALPS::utilities` owns utility symbols and links Boost.Filesystem and platform threads without the simulation runtime, HDF5 or BLAS/LAPACK.

`ALPS::hdf5` owns archive symbols, exception exports, shared archive state and the NGS signal handler that closes archives. It links utilities, HDF5, Boost.Filesystem, Boost.Thread and platform threads. Archive and signal code remain together to preserve cleanup behavior. A parallel HDF5 provider can bring its own MPI dependency.

`ALPS::params` owns typed values, proxies, checkpoint I/O and the `paramvalue_source` interface used by Python bindings. It links HDF5 and Boost.Serialization; MPI builds additionally link MPI and Boost.MPI. Its parameter-file constructor, XML input and older `Parameters` conversion remain in `params/adapters/`, compiled into `ALPS::alps`. The two conversion headers now live in `params/adapters/include/alps/ngs/`, retaining their public names through the aggregate interface rather than the typed params target. The older `alps::Parameters` implementation itself belongs to `legacy_parameters/`.

`ALPS::osiris` owns dump/process APIs, XDR symbols and the communication state used by `comm_init()` and `is_master()`. It links Boost.Serialization/Filesystem and, when enabled, MPI. This gives communication state one owner and permits MaxEnt to preserve existing diagnostic gating without linking the simulation runtime.

`ALPS::xml` owns XML parsing, attributes, handlers, output streams and stylesheet lookup, with Boost.Filesystem and Boost.Regex dependencies. XML file-to-parameter conversion still belongs to the params adapters in `ALPS::alps`. The separate `plotting/` header owner keeps `<alps/plot.h>` and its older `Parameters` dependency outside the XML component.

`ALPS::cli` owns the existing `mcoptions` and `parseargs` implementations, linking utilities and Boost.ProgramOptions. Their installed headers remain `<alps/ngs/mcoptions.hpp>` and `<alps/parseargs.hpp>`. The two existing option grammars, defaults, filename rules and error behavior are preserved; these parsers do not read parameter files.

`ALPS::maxent` links params, HDF5, utilities, Osiris and its Boost/numerical providers. It owns its deterministic run loop; its numerical calculations and stop-callback ordering are preserved. The `maxent` executable adds `ALPS::cli` for the existing `mcoptions` parser and reads typed params from HDF5 directly. Neither the solver nor executable links `ALPS::alps`. Its public callable API remains `<alps/solvers.hpp>`.

`ALPS::alps` links the extracted runtime components publicly. These libraries follow `BUILD_SHARED_LIBS`, with component-specific generated export headers. Python extensions require shared runtime libraries and package one copy of every component in `ALPS_RUNTIME_TARGETS`. Rebuild downstream binaries after the XML and CLI extractions, as after the earlier library splits.

Physical ownership does not imply independent linkability for every module. `ALPS::headers` still exposes the aggregate compile interface for the simulation modules, including the remaining expression/older-parameters cycle. The foundation interfaces have narrower dependencies, but package discovery still checks the full SDK dependency set. Shared build policy remains in the root CMake files and `cmake/`.

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

With application builds enabled, `alps-module-architecture.json` inventories 35 source owners, 517 public include spellings and 654 production files. The owners include `numeric_io`, `cli`, `plotting`, separate MaxEnt solver/executable owners, and eight [tool groups](../tools/README.md). Tool ownership includes historical inactive C++ sources without adding executable targets. These are ownership counts, not counts of independent libraries or passing tests. The earlier code checkpoint `f6f4501c0` had 24 owners, before the CLI and plotting modules were separated.

The foundation include cycle involving containers, HDF5, numerics, utilities and XML is removed. The current observed include graph retains the separate two-module cycle between `expression` and `legacy_parameters`; it still needs deliberate reconciliation. Dependency declarations constrain new include edges. Regenerate the report after changing module ownership or dependencies.

The report inventories 82 exact unresolved file/include pairs already present at baseline `0f7b995d5`: 80 in dormant accumulator code and two in the optional `USE_LATTICE_CONSTANT_2D` graph backend. Each exemption names its file, include and reason; they do not establish support for those inactive paths. Resolve or remove these dependencies deliberately rather than adding broad exclusions.

With `ALPS_BUILD_TESTING=ON`, CTest runs `module_architecture` and writes `<build-dir>/alps-module-architecture.json`; this requires a Python interpreter ≥ 3.10. Builds with testing disabled do not need Python for module configuration or manifest generation.

## Validation

`ALPS_BUILD_TESTING` controls module-local and central native tests. `ALPS_BUILD_APPLICATIONS` additionally controls MaxEnt and its regression. After building, focused existing labels include:

```sh
ctest --test-dir <build-dir> --output-on-failure -L '^(utility|hdf5|params|osiris|parser|cli|maxent)$'
```

Run the full native suite for shared-interface changes. The SDK consumers in `tests/cmake/consumer/` exercise aggregate and component links; integration tests for older inputs and observables keep `ALPS::alps`. Check shared/static consumers, installed-header ownership, SDK relocation and Python extension interoperability. The CLI contract checks existing defaults, option spellings, filename rules, help and execution-mode handling. The MaxEnt regression also exercises early stop, callback exceptions and completion ordering. These are validation requirements; this ownership map does not assert that every platform/configuration has passed them.
