# Proposal: TOML configuration, shared validation, and parameter sweeps

**Status:** Draft for discussion; no implementation or adoption is implied.

**Assessment baseline:** ALPS `master` at
[`1932d6385`](https://github.com/ALPSim/ALPS/commit/1932d63854d6e5d1a0f44c19a2b2bede45400b86).

## Overview

I propose a new configuration and experiment architecture for ALPS, with incremental
backend migration. It would use TOML input, one typed configuration model, shared
schemas, and validation at defined preparation and initialization stages. Parameter
sweeps and sampling replicas would be first-class features. HDF5 would remain the
numerical result format and the target for migrated checkpoint implementations.

The main change is to separate responsibilities that currently overlap:

| Responsibility | Proposed owner |
| --- | --- |
| Parse authored configuration | TOML reader or legacy input adapter |
| Describe types, defaults, and constraints | Shared framework and solver schemas |
| Define runs, sweep coordinates, replicas, and assigned seeds | Experiment planner and persistent experiment definition |
| Track execution progress, scheduling, and retries | Existing execution backend, accessed through an adapter |
| Resolve expressions and validate scientific settings | Configuration resolver and staged solver validators |
| Hold effective settings for one run | Owned, typed `RunConfig` |
| Hold counters, RNG state, and measurements | Mutable solver state |
| Store numerical data and restart state | Separate result and checkpoint writers |

I recommend `toml++` for native TOML parsing and serialization. Its types should
remain internal to the reader/writer. The initial design would reuse the existing
expression engine and HDF5 archive implementation through explicit boundaries.

The experiment model owns what should run; the execution backend owns its progress.
The new configuration model is the long-term authority. Both existing parameter
containers become compatibility boundaries as solvers migrate. A format converter
can support that transition, but the first milestone must also migrate one solver's
configuration consumption to the new typed interface.

That milestone should cover one solver, provisionally `spinmc`, including single
runs, product and zip sweeps, staged validation, effective configuration output,
and integration with its existing execution and restart mechanisms. Existing input
formats and Python APIs remain operational during migration. Uniform checkpoint
storage across all backends and content-derived run identities are deferred.

The decisions needed before implementation are the schema ownership model,
expression and sweep semantics, compatibility commitments, and first solver.
All APIs, commands, and new-format examples below are illustrative proposals.

## Contents

- [Current behavior and motivation](#current-behavior-and-motivation)
- [Goals and scope](#goals-and-scope)
- [Configuration architecture](#configuration-architecture)
- [TOML format and dependencies](#toml-format-and-dependencies)
- [Schemas and validation](#schemas-and-validation)
- [Defaults, overrides, and paths](#defaults-overrides-and-paths)
- [Scientific expressions](#scientific-expressions)
- [Parameter sweeps and experiment planning](#parameter-sweeps-and-experiment-planning)
- [Results, checkpoints, and provenance](#results-checkpoints-and-provenance)
- [Python and MPI interfaces](#python-and-mpi-interfaces)
- [Migration and first milestone](#migration-and-first-milestone)
- [Decisions for discussion](#decisions-for-discussion)
- [References and provenance](#references-and-provenance)
- [Short issue text](#short-issue-text)

## Current behavior and motivation

ALPS currently distributes configuration responsibilities across several APIs:

- [`alps::Parameters`](../../src/alps/parameter/parameters.h) stores values as
  strings and supports the legacy parameter language and XML.
- [`alps::params`](../../src/alps/ngs/params.hpp) introduces typed values, but its
  [text-file constructor](../../src/alps/ngs/lib/params.cpp) still uses the legacy
  parser. Language bindings add further value and decoding behavior.
- Defaults, conversions, and scientific constraints appear in solver
  constructors, for example in
  [`AbstractSpinSim`](../../src/apps/mc/spins/abstractspinsim.h).
- [`Parameters::save()`](../../src/alps/parameter/parameters.C) evaluates
  expressions and chooses numeric representations while writing HDF5.
- [`ParameterListParser`](../../src/alps/parameter/parameterlist_p.h) supports
  evolving global defaults and local parameter blocks. A migration must preserve
  the resulting configurations, not merely translate individual assignments.
- Python [`writeInputFiles()`](../../python/pyalps/src/pyalps/tools.py) combines
  seed assignment, mutation of parameter dictionaries, job generation, and XML
  writing.
- [Python result readers](../../python/pyalps/src/pyalps/load.py) and
  [native restart code](../../src/alps/mcbase.cpp) rely on existing HDF5 paths and
  encodings.
- [Parapack tasks](../../src/alps/parapack/job.C) already own clone progress,
  continuation, and aggregation. Its
  [seed initialization](../../src/alps/parapack/clone_info.C) distinguishes worker
  and disorder seeds and can depend on process topology.
- Persistence differs between backends. For example,
  [Parapack clone persistence](../../src/alps/parapack/clone.C) can write HDF5
  metadata and measurements alongside separate XDR worker state.

This makes supported parameters and their defaults difficult to discover. It
also complicates validation before execution, consistent diagnostics across
interfaces, and reconstruction of the effective settings behind a result.

The existing system already has important compatibility behavior. In particular,
[Python checkpoint tests](../../tests/pyalps/test_checkpoint_contracts.py) require
mutable parameter objects, arrays, nested metadata, and exact large integers.
Those requirements must inform the migration rather than be discarded by a new
primitive-only parameter container.

## Goals and scope

I propose the following goals:

1. A documented, human-readable input format.
2. One authoritative definition of parameter types, defaults, and constraints.
3. Consistent configuration semantics across native and Python entry points.
4. Early actionable input checks, followed by model-aware initialization checks.
5. Explicit expressions, sweep expansion, replica identity, and seed assignment.
6. A reproducible record of effective settings and resource identities.
7. Versioned result and restart contracts with a gradual compatibility path.

The replacement architecture should be established from the first implementation;
adapters allow its adoption one solver and backend at a time. The scope is divided
as follows:

| Commit to in the initial design | Migrate or extend later |
| --- | --- |
| One typed configuration authority and shared validation rules | Remaining solver constructors and legacy parameter APIs |
| Structured TOML with explicit legacy-name and expression mappings | General XML model/lattice API migration |
| Product/zip sweeps, persisted run IDs, and recorded seed assignments | Advanced sweep composition and plan revision tooling |
| Immutable prepared settings and explicit resume-time changes | Complete provenance dependency graphs and optional configuration hashes |
| Versioned result/checkpoint metadata and an artifact-aware loader | Uniform checkpoint storage, publication, and migration across backends |
| A solver consuming typed settings through an existing execution backend | Additional backend adapters and generated editor schemas |

Scientific algorithms, a replacement cluster scheduler, and an HDF5-wrapper
rewrite are outside this proposal's initial implementation. Compatibility adapters
must have explicit ownership and a migration path; the new model should not become
a third permanent parameter system beside the existing two.

## Configuration architecture

```text
TOML file / Python mapping / legacy input adapter
                         |
         Parsed document with source information
                         |
         Apply overrides and expand experiment
                         |
         Resolve preparation-time defaults and expressions
                         |
         Check inputs and assign run IDs, seeds, and outputs
                         |
         Persistent experiment definition of prepared runs
                         |
         Existing execution backend through an adapter
                         |
         Initialize model/lattice context; finish validation
                         |
         Solver using final typed RunConfig and mutable state
                         |
          Separate result and checkpoint writers
```

The configuration core should be usable without Python, MPI initialization, or
an HDF5 archive. TOML parsing, legacy conversion, persistence, and language
bindings should be separate adapters to that core.

The parsed document should retain each value's origin: input file and position,
Python mapping path, command-line argument, or schema default. The initial
provenance record includes effective values, applied defaults and overrides, and
resource identities. A full dependency graph for derived values is optional later.

Prepared configuration should be immutable through its public API. Context-dependent
defaults and expressions remain explicitly represented until initialization can
resolve them. Initialization produces a final typed `RunConfig` and the required
scientific context; it does not silently rewrite the authored experiment definition.
The effective snapshot records these initialization decisions.

Migrated solvers consume typed fields. Dynamic model parameters can remain a
validated typed map or symbolic representation owned by the model component.
Counters, RNG state, measurements, and adaptive algorithm state belong elsewhere.
Adapters may translate settings into legacy containers for dependencies that have
not migrated, but migrated fields must have only one owner for defaults and checks.

Parsing and basic validation must not open simulation outputs or construct a
solver. Model-aware checks may need graph construction or other scientific
initialization. Those checks must complete before simulation steps begin, and
runtime invariants remain the solver's responsibility. A validation report must
state which stages ran and which checks require execution-time initialization.

An illustrative API boundary is:

```cpp
auto document = alps::config::read("experiment.toml");
auto prepared = alps::config::prepare(document, overrides, solver_registry);
if (!prepared)
    return report_errors(prepared.errors());
execute(prepared.value(), execution_options);
```

These names describe responsibilities, not a commitment to a particular C++
error-handling or container implementation.

## TOML format and dependencies

I recommend [`toml++`](https://github.com/marzer/tomlplusplus) for native reading
and writing. It fits ALPS's C++17 baseline, has no required runtime dependencies,
and provides source locations for diagnostics. The parser should remain an
implementation dependency; public solver headers should not expose `toml::table`.

The first format version should use a documented TOML 1.0 profile with experimental
language extensions disabled. Dependency versions should be pinned and available
to offline builds, with supported system-package integration. Bundled dependency
notices must be retained.

[`toml11`](https://github.com/ToruNiina/toml11) is an alternative if preserving
comments while editing authored files is a core requirement. The initial workflow
would leave authored files intact and write separate resolved snapshots, which
reduces that requirement.

An illustrative single-run input is:

```toml
format_version = 1
solver = "spinmc"

[lattice]
name = "chain lattice"

[lattice.parameters]
L = 64

[model]
name = "Heisenberg"

[model.parameters]
J = 1.0

[monte_carlo]
temperature = 0.5
sweeps = 100_000
thermalization = 10_000
seed = 42

[output]
results = "results.h5"
checkpoint = "checkpoint.h5"
```

Framework settings have explicit namespaces. Model/lattice symbols such as `J`
and `L` remain within the corresponding parameter tables. Existing XML libraries
remain usable through references; the chosen solver adapter determines how those
references map to current interfaces.

Each migrated solver must document its legacy mappings, for example
`monte_carlo.temperature` to `T`, `monte_carlo.sweeps` to `SWEEPS`, and
`lattice.parameters.L` to `L`. These mappings also determine compatibility fields
in result metadata used by existing analysis scripts. Ambiguous mappings must fail
explicitly; nested tables must not be flattened by discarding their namespaces.
New code uses the structured names while compatibility interfaces retain the
scientific vocabulary expected by existing users.

Unknown format versions must fail clearly. TOML dates and times should be rejected
unless explicitly admitted by a field schema. Complex values require a documented
schema representation, such as `{ real = 1.0, imag = 0.5 }`; they must not be
inferred from arbitrary strings. Large numerical arrays should normally be
referenced through typed datasets rather than embedded in configuration files.

The TOML profile must specify the signed 64-bit integer boundary. Values outside
that range in historical archives or Python objects need an explicit compatibility
representation or a checked rejection when entering the new API. No path may
silently narrow them or convert them through floating-point storage.

## Schemas and validation

Each solver should provide a schema composed from shared framework fields and
solver-specific fields. Each field definition should declare:

- Canonical path and accepted type.
- Whether the field is required.
- Literal or computed default.
- Validation/resolution stage and component responsible for the field.
- Allowed values, bounds, and shape constraints.
- Scientific meaning and unit convention.
- Deprecated aliases and compatibility mappings.
- Whether it may change when resuming a checkpoint.

I recommend C++ schema descriptors associated with typed fields, plus named
functions for scientific constraints. The first implementation should establish
field ownership, checked decoding, defaults, and validation. Generating reference
documentation and editor-facing JSON Schema can follow once those contracts work
for real solvers. A general-purpose reflection framework is unnecessary.

Python should invoke the native validation implementation rather than maintain
independent scientific rules. External JSON Schema is an alternative authority
for structural rules, but it would still require scientific validators,
expression resolution, and source-location mapping. If chosen, its supported
schema dialect and default-application behavior must be explicit.

Validation must distinguish missing values from malformed values, booleans from
integers, integers from fractional values, and checked conversion from overflow
or precision loss. Defaults apply only to absent fields. Finite-value rules should
be field-specific; infinity may be meaningful for some scientific settings.

For example, `sweeps = 1000.5` should fail rather than truncate. A real-valued
field may accept an integer only under a documented conversion policy. Unknown
framework or solver fields should fail with suggestions where possible.

Validation should have three explicit stages:

1. **Input checks:** syntax, field names, required values, checked types, ranges,
   and constraints whose inputs are available during preparation.
2. **Initialization checks:** constraints and defaults requiring the selected
   graph, model, or algorithm representation. They finish before simulation steps.
3. **Runtime invariants:** properties that depend on evolving simulation state.

For example, [`spinmc` chooses its default update method from graph frustration
and coupling properties](../../src/apps/mc/spins/abstractspinsim.h). Resolving this
default requires scientific context. Moving the rule to a named initialization
function shared by the solver and an optional thorough validation mode avoids
duplicating it in an incomplete generic preflight validator.

Other rules include temperature versus inverse temperature and compatible
field/update combinations. A generic rule such as requiring thermalization to
be less than measurement sweeps should not be imposed without establishing that
solver's semantics.

Model/lattice XML defaults are not complete schemas. Symbols can be introduced by
expressions or site/bond-type substitutions such as `J#`. Model-aware validators
must account for these mechanisms. Strict unknown-field rejection applies where
a component declares a complete schema. A legacy model adapter must state its
validation coverage, preserve its established symbol resolution, and report
unresolved symbols at the appropriate stage. It must not silently claim full
validation or use extension data to bypass checks of known solver controls.

Independent failures should be collected where possible:

```text
run.toml:18:10: monte_carlo.sweeps
  expected a positive integer; received 1000.5

run.toml:19:1: unknown field "thermalisation"
  did you mean "thermalization"?
```

Diagnostics should also expose stable error codes, configuration paths, and run
coordinates for Python and automation. Errors that prevent further evaluation
should suppress dependent diagnostics rather than produce misleading cascades.

## Defaults, overrides, and paths

For a single run or explicit run list, precedence would be:

```text
schema defaults -> base input -> per-run overrides -> explicit CLI overrides
```

Computed defaults and expressions must use final effective inputs at their declared
resolution stage. Canonical-name and alias conflicts within the same layer should
fail. Intentional replacement between layers should retain the winning origin and
the applied override.

Sweep axes are an explicit additional layer described below. CLI overrides of an
axis target should be rejected rather than silently disabling the sweep.

There should be no implicit environment substitution in the new format. If added,
environment-derived values should be explicitly requested and their effective
values recorded. Legacy adapters retain the existing substitution semantics.

File-supplied paths should resolve relative to their source file; CLI-supplied
paths should resolve relative to the invocation directory. Execution should not
depend on a worker's incidental working directory. A portable plan should retain
resource identity and allow explicit relocation with identity checks.

## Scientific expressions

Scientific expressions are an existing capability that the migration must
preserve. The new format should distinguish them from literal strings:

```toml
[constants]
J0 = 1.0

[model.parameters]
J = { expr = "2 * constants.J0" }
```

Configuration expressions should use a defined namespace, detect dependency
cycles and undefined symbols, and resolve during preparation when their inputs
are known. Context-dependent expressions resolve during initialization or in their
declared scientific context. Errors should point to the authored expression. Both
the expression and effective value should be retained when applicable.

Site-, bond-, and quantum-number-dependent expressions must remain validated
symbolic objects until their evaluation context exists. Known context symbols
must be distinguished from misspelled configuration symbols.

I recommend reusing the existing expression engine through an adapter initially.
Its current [identifier grammar](../../src/alps/parser/parser.h) does not accept
dotted configuration paths. Supporting `constants.J0` therefore requires an explicit
reference resolver or parser integration, tested alongside legacy scientific
symbols. The expression syntax is part of the new configuration API, with a
documented mapping to the existing evaluator; a TOML parser alone does not supply it.
Random expressions and disorder generation need explicit RNG ownership and
evaluation timing. Saving a configuration must never evaluate expressions again
or consume random numbers.

## Parameter sweeps and experiment planning

Here, a **parameter sweep** means executing multiple configurations. It is
distinct from `monte_carlo.sweeps`, which controls work within one simulation.

I propose supporting individual runs, explicit run lists, and parameter sweeps
with sampling replicas. Planning should produce a persistent experiment definition
with prepared runs and the status of deferred initialization checks. This is the
public experiment interface; backend-specific XML jobs or other launch artifacts
may be generated by adapters during migration.

### Product and zip expansion

An illustrative product experiment is:

```toml
format_version = 1
solver = "spinmc"

[lattice]
name = "chain lattice"

[lattice.parameters]
L = 32

[model]
name = "Heisenberg"

[model.parameters]
J = 1.0

[monte_carlo]
sweeps = 100_000
thermalization = 10_000

[experiment]
name = "heisenberg-temperature-scan"
base_seed = 42
sampling_replicas = 4

[sweep]
mode = "product"

[[sweep.axes]]
parameter = "lattice.parameters.L"
values = [32, 64, 128]

[[sweep.axes]]
parameter = "monte_carlo.temperature"
values = [0.25, 0.5, 1.0, 2.0]

[output]
directory = "results"
```

This defines 12 parameter combinations and four sampling replicas each: **48 runs**.

Two modes should be available initially:

- `product`: Cartesian product of all axis values.
- `zip`: Combine corresponding entries; every axis must have the same length.

For example, this fragment defines three configurations:

```toml
[sweep]
mode = "zip"

[[sweep.axes]]
parameter = "lattice.parameters.L"
values = [32, 64, 128]

[[sweep.axes]]
parameter = "monte_carlo.sweeps"
values = [100_000, 200_000, 400_000]
```

There should be no inferred expansion mode or implicit broadcasting. The format
must specify deterministic enumeration order for initial preparation. Persisted
IDs and assigned seeds survive subsequent scheduling reordering; recreating a
plan after editing or reordering its axes need not preserve those assignments.

Ordinary arrays elsewhere remain parameter values. For a vector-valued axis,
each entry of `values` represents a complete vector candidate, checked against
the field's type and shape. Empty axes, duplicate targets, invalid candidate
types, and unknown targets should fail clearly. An axis may supply a required
field absent from the base input.

Axis paths refer to configuration structure, not legacy flat parameter names.
The path grammar must support quoted segments for model symbols containing dots
or other special characters; a plain string split on dots is insufficient.

The first implementation should allow declared scientific and solver fields as
targets. Solver selection, format versions, replica counts, seeds, and output
destinations should not themselves be axes. Combining an explicit run list and
a sweep in one input should fail until composition semantics are specified.
Grouped product/zip composition and range shorthand can follow later.

### Resolution and invalid combinations

The preparation order should be:

1. Load base settings and apply experiment-wide overrides.
2. Expand axis combinations, replacing base values at their target paths.
3. Resolve preparation-time defaults and expressions for each combination.
4. Run input checks and record any required initialization checks.
5. Expand sampling replicas and assign run IDs, seeds, and destinations.
6. Check the complete definition for conflicts and publish it only if preparation
   succeeds. Initialization checks still gate simulation startup for each run.

Expressions therefore see each combination's axis values. CLI overrides targeting
an axis must fail with an explanation. Non-axis overrides are applied before
derived values are resolved. The base input may omit fields supplied by axes;
required-field validation must account for this.

Invalid combinations must not be silently skipped. Preparation failures should
identify the coordinates and constraints involved. A combination that fails a
deferred initialization check must be reported as a failed run through the backend;
it must not disappear from the experiment. An eventual filtering feature must be
explicit and record excluded combinations with reasons.

### Identity and seeds

Initial preparation should assign and persist an experiment ID, parameter-point
IDs, and run IDs. Each run also records its sweep coordinates and sampling replica
index. Explicit user IDs or simple generated IDs are sufficient. Once a definition
is saved, execution reordering must preserve its IDs and assigned seeds.

Configuration hashes may be added as provenance later. They should not determine
run identity or seeds in the first implementation. Changing a permitted sampling
budget on resume keeps the run ID and records the validated setting change.
Simple duplicate parameter points should produce warnings; deliberate repeated
points remain legal and receive distinct run IDs. Semantic equivalence analysis
of arbitrary expressions is unnecessary.

The initial seed policy should assign task-level seeds once from the recorded base
seed and deterministic preparation order, using a documented policy identifier.
Explicit run lists may provide task seeds directly. Mapping into each backend's
supported seed range must be checked. Persist the actual assignments so execution
and restart do not recompute them. Conflicting seed declarations must fail clearly.

The architecture must distinguish these scientific identities:

| Identity | Meaning |
| --- | --- |
| Parameter point | One combination of scientific and sampling settings |
| Disorder realization | A particular random physical instance, when applicable |
| Sampling replica | A separate sampling chain for the selected physical instance |
| Worker/rank | An execution participant within a run |

Sampling replicas must not implicitly request new disorder realizations. A disorder
capable adapter must expose an explicit realization/seed policy, including whether
the same realization is shared across parameter points. Replica statistics and
disorder averages must retain their distinct identities. The initial clean-model
`spinmc` prototype covers sampling replicas; unsupported disorder requests must
fail explicitly rather than be ignored.

Each backend adapter must document the mapping to its tasks/clones and its worker
and disorder seeds. It must account for backend-created chains so that four planned
sampling replicas do not silently expand into four times the backend's automatic
clone count. Unsupported execution layouts should fail before submission.
Existing Parapack worker streams can depend on process topology.
The first guarantee is stable persisted task-level seed assignments under execution
reordering. Internal-stream invariance across process counts or backend changes is
not assumed, and distinct seeds alone do not prove statistical independence.

Automatic base-seed generation, if supported, happens once during preparation and
is persisted. Legacy imports preserve resolved historical seeds and record the
backend seed policy. They must not silently adopt a new derivation scheme.

### Inspection and persisted plans

Illustrative commands are:

```sh
alps plan experiment.toml
alps plan experiment.toml --write plan.toml
alps run --plan plan.toml
```

Planning should report axes and candidate counts, combinations, replicas, total
runs, prepared settings, seeds, destinations, and validation errors. It must also
identify automatic values and checks awaiting scientific initialization. A preview
must not create simulation output files. An optional thorough validation mode can
initialize the required scientific context without taking simulation steps.

The planner should calculate expansion size using overflow-checked arithmetic
before materializing combinations and enforce a documented configurable run-count
limit. Explicit value lists are sufficient initially; future range generators
must define endpoints and floating-point stepping precisely.

A persisted definition contains prepared run settings, explicit unresolved
context-dependent values, resource identities, schema/resolver versions, and seed
assignments. Execution uses it without re-expanding edited source inputs. Deferred
decisions must use a compatible resolver; changed default semantics require an
explicit migration or rejection. Final effective settings are recorded after
initialization and become authoritative for restart. Definition loading validates
compatibility and any relocated resources.

### Outputs, execution, and experiment restart

Each run needs unambiguous artifact locations. One possible layout is:

```text
results/
  plan.toml
  runs/
    run-000001/
      effective.toml
      results.h5
      checkpoint.h5
    run-000002/
      ...
```

The directory layout is optional and should not define identity. An adapter may
retain backend filenames while recording their locations and roles in the
experiment definition or associated backend artifact index. Existing destinations
require an explicit resume or replacement policy.

Each result should store experiment/parameter-point/run IDs, sampling replica and
disorder identities, coordinates, and effective settings. A new experiment loader
must enumerate declared result artifacts and check their roles. The existing
[`getResultFiles()`](../../python/pyalps/src/pyalps/tools.py) can fall back to matching
all HDF5 files, which would include checkpoints in the example above. Legacy
discovery remains available for old workflows; new-format discovery must not rely
on that heuristic. Existing per-file result readers should still be usable where
their dataset contract is preserved.

Replica identities should remain available in stored data even if a backend
publishes an aggregate. Existing clone evaluators can be reused through an explicit
mapping that records which chains and disorder instances contributed. Aggregation
is a scientific operation and must not conflate sampling and disorder averages.

Scheduling backends should consume the same prepared definitions. They must
distinguish concurrency between runs from parallelism within a run: an MPI solver
may need several ranks for one configuration.

The existing backend owns progress, completion, and retry state. The adapter maps
persisted run IDs to backend task/clone identities and exposes status through a
common view. That view must not become a second independently writable progress
database. Submission/resume must reconcile the saved mapping to avoid launching
the same prepared run twice.

Resuming an experiment validates its definition and settings, then delegates
continuation and completion decisions to the selected backend. Backend recovery
limits must be documented and tested. Permitted changes such as extending a
sampling budget retain the run ID and seed assignment and produce a recorded
configuration update. Sophisticated plan revisions and cross-backend checkpoint
migration can follow later.

## Results, checkpoints, and provenance

| Artifact | Representation | Contract |
| --- | --- | --- |
| Authored configuration | TOML | User-maintained settings, expressions, experiment definitions |
| Prepared experiment | Persisted typed definition, with TOML export | Run IDs, coordinates, prepared settings, assigned seeds, deferred checks |
| Effective configuration | Generated TOML and embedded HDF5 metadata | Final initialization decisions, effective values, defaults, overrides, resource identities |
| Results | HDF5 | Typed numerical outputs and analysis metadata |
| Checkpoints | Backend-specific during migration; HDF5 for migrated implementations | Complete restart state with explicit backend and format identity |

Serialization should preserve typed values and explicit deferred expressions
without evaluating them or inferring their types. The experiment definition owns
prepared intent; the final effective configuration owns the settings actually used
by a run. Their relationship and versions are recorded. Human-readable exports
are snapshots, and conflicting copies must not be independently accepted as
authorities for restart.

Results should record ALPS/solver and schema versions, model/lattice resource
identities, effective configuration, RNG information, and experiment coordinates.
New artifacts need explicit roles such as result or checkpoint, plus format and
backend identifiers. Result and checkpoint formats need separate version contracts:
being able to read a result does not imply being able to resume it. Resource hashes
are useful provenance where available, but complete semantic configuration hashing
is not required for initial adoption.

The checkpoint contract requires full RNG, accumulator, and solver state, which
may occupy multiple backend artifacts during migration. The adapter must identify
the complete set; an HDF5 measurements file alone must not be advertised as a
resumable checkpoint. Restart compares effective settings with requested changes
using shared field policies and the backend's actual capabilities. A field declared
mutable is permitted only if that backend can safely apply the change. New defaults
must not reinterpret saved state, and unsupported versions must fail clearly.

The first backend migration should preserve its existing restart behavior and
document publication/recovery limitations. As checkpoint implementations migrate,
they should publish closed, complete generations, retain a recoverable prior
generation, and validate into temporary state before replacing live state. Atomic
visibility and durability must be specified separately. Establishing that target
contract now does not require rewriting every backend's storage in the first
milestone, or treating existing multi-file writes as already transactional.

Initially preserve established HDF5 numerical paths and encodings, including
`/simulation/results`, and provide legacy parameter-name projections where analysis
needs `/parameters/T` or similar fields. Structured configuration has one authority;
compatibility projections are generated and checked. New metadata must not assume
that old discovery utilities understand artifact roles. The new loader handles
classification while compatible numerical readers can be reused. Independent
readers should be able to inspect supported numerical datasets without
ALPS-specific decoding.

A future archive-wrapper replacement, for example evaluating
[HighFive](https://github.com/highfive-devs/highfive), should be a separate change
with format-compatibility tests rather than a prerequisite for TOML support.

## Python and MPI interfaces

New Python preparation APIs should take mappings or file inputs through the same
normalization and validation logic as native executables. Preparing a run should
create an owned snapshot; later dictionary or array mutations should not alter it.

Existing mutable `ngs.params` semantics remain available through compatibility
APIs. Arrays, nested metadata, and exact integers need documented mappings at the
new boundary. Unsupported values must fail explicitly rather than be stringified
or narrowed. Metadata should be separated from solver settings so arbitrary
annotations do not bypass scientific validation.

Python TOML utilities can use `tomllib` on Python 3.11+, Tomli on older supported
versions, and Tomli-W for generation. TOML Kit is an option for formatting-preserving
editing. File inputs requiring identical native diagnostics should use the native
reader, and any alternative parser path must pass the same format-profile tests.

Python should expose the shared configuration API, experiment planning, and loading
of declared result artifacts with their coordinates, sampling replicas, and disorder
identities. Input preparation should work without launching simulations and report
deferred scientific checks. It must not introduce a separate Python defaults or
validation implementation.

For MPI, the designated rank should prepare settings and distribute the typed
representation. Initialization checks must be coordinated across the participating
ranks before simulation steps. Failure must terminate startup coherently rather
than leave peers blocked in collectives. Existing execution backends retain
communicator and scheduling ownership, distinguishing parallelism within runs from
distribution between runs. HDF5 writes need a designated writer or explicitly
supported parallel-I/O path. Adapters must record topology-dependent seed behavior
rather than promise internal-stream invariance they cannot provide.

## Migration and first milestone

I propose the following sequence:

1. Establish fixtures for legacy parsing, expressions, evolving parameter-list
   defaults, seeds, numeric precision, Python behavior, and checkpoints.
2. Define the authoritative typed configuration model, schemas, preparation and
   initialization validation interfaces, and TOML reader. Put legacy conversions
   at explicit compatibility boundaries.
3. Migrate one solver's configuration consumption, provisionally `spinmc`, to typed
   fields. Move ownership of those defaults/checks out of the legacy path and reuse
   the shared definitions for compatibility callers.
4. Add product/zip preparation, persistent experiment definitions, recorded seeds,
   and an adapter to that solver's existing execution backend. Add result-role
   metadata, effective configuration output, and artifact-aware loading.
5. Expose shared preparation through Python and add legacy conversion reports.
   Migrate further solvers and backends with a compatibility/capability matrix.
6. Migrate checkpoint implementations to the defined persistence contract
   incrementally. Agree on API deprecations once compatibility coverage is adequate.

This is a vertical implementation of the new architecture. Translating TOML into
legacy parameter text may help integration but does not by itself complete the
milestone. Scientific dependencies still using old parameter containers can remain
behind adapters while the chosen solver consumes typed configuration directly.

`spinmc` is a useful scheduler integration target, but its
[model factory](../../src/apps/mc/spins/spinmc_factory.C) dispatches classical models
directly. Its success will not validate every XML Hamiltonian or disorder workflow.
Include a focused general-model expression/parameter adapter fixture and schedule
a general-model solver as a subsequent migration target before claiming broad
schema coverage.

Legacy conversion should materialize explicit resolved run definitions when the
original syntax cannot be represented faithfully as a sweep. It must report
ambiguous or unsupported transformations. TOML support must not change the meaning
of existing text files or infer the format from whether one parser happens to fail.

The first milestone should satisfy these acceptance criteria:

- A representative single-run TOML input executes a solver consuming the new typed
  configuration interface; migrated defaults and checks have one authoritative owner.
- Input-only planning creates no simulation outputs and reports deferred checks.
  Initialization checks gate simulation startup, and failures are actionable.
- The saved experiment records IDs, prepared settings, and task seeds. The effective
  snapshot records final initialization decisions and the settings actually used.
- Representative legacy and TOML inputs produce equivalent effective parameters
  and scientific results under the same execution conditions.
- Product and zip sweeps expand correctly; malformed axes and conflicts fail.
- Expressions and computed defaults resolve per combination at their declared stages;
  namespaced references and legacy scientific symbols have tested mappings.
- Vector-valued parameters remain unambiguous, and numeric boundary cases are checked.
- Persisted run IDs and task seed assignments survive scheduling reordering and
  allowed resume-time budget changes. Duplicate points receive distinct IDs.
- Sampling replicas have an explicit backend mapping. Unsupported disorder or
  topology requirements fail clearly; no stronger RNG guarantee is implied.
- Run destinations are unique. Results identify their role, coordinates, replicas,
  and any disorder realization. New discovery excludes checkpoint artifacts.
- A small experiment executes and resumes through the chosen backend without a
  second progress database, duplicate submission, or changed seed assignments.
- A deterministic restart fixture preserves scientific behavior under the chosen
  backend's supported conditions. Incompatible settings are rejected before state
  loading; storage/recovery limitations are recorded in the capability matrix.
- Existing per-file result readers and Python compatibility tests continue to pass;
  new experiments are loaded through the artifact-aware entry point.
- Installed SDK consumers do not need TOML implementation types in solver headers.

Existing parameter, native/Python checkpoint, and HDF5 regression suites should be
reused. New tests should focus on public contracts and scientific equivalence.
Parallel adapters need explicit multi-rank failure and restart tests as they
migrate. Generation-retaining checkpoint publication, transactional state loading,
full provenance graphs, generated editor schemas, and advanced plan revision tools
are follow-up milestones rather than prerequisites for the first solver path.

## Decisions for discussion

1. Is the ownership split between the new configuration/experiment model and
   existing execution backends appropriate?
2. Is `toml++` preferred, or is comment-preserving editing important enough to favor
   `toml11`?
3. Should native schema descriptors be authoritative, and how should staged
   validation and third-party solver/model extensions be registered?
4. Are the structured names, legacy mappings, expressions, paths, and product/zip
   semantics suitable?
5. Is the separation of parameter points, disorder realizations, sampling replicas,
   and workers sufficient? Which backend task-seed policy should the first adapter use?
6. Which historical inputs, Python APIs, and checkpoint versions must remain
   supported, and for how long?
7. Is `spinmc` the best first target, and which general-model solver and persistence
   backend should migrate next?

## References and provenance

This is an architectural proposal based on the linked ALPS source at the stated
baseline. Examples illustrate proposed ALPS interfaces; no external implementation
is copied into this document and no dependency is added by this branch.

- [TOML 1.0 specification](https://toml.io/en/v1.0.0)
- [toml++](https://github.com/marzer/tomlplusplus) and its
  [source-region API](https://marzer.github.io/tomlplusplus/structtoml_1_1source__region.html)
- [toml11](https://github.com/ToruNiina/toml11)
- [Python TOML documentation](https://docs.python.org/3.11/library/tomllib.html)
- [nlohmann-based JSON Schema validator](https://github.com/pboettch/json-schema-validator),
  an option if external structural schemas are selected
- [HighFive](https://github.com/highfive-devs/highfive), a possible separate archive
  modernization topic

Implementation should preserve provenance and license notices and add relevant
scientific credit through `CITATION.cff` and `CITATIONS.yaml` following the
[contribution policy](../../CONTRIBUTING.md#provenance-and-scientific-credit).
This discussion document does not change scientific implementations or citation
ownership.

## Short issue text

**Suggested title:** Proposal: Modernize ALPS configuration, validation, and parameter sweeps

ALPS currently spreads parameter parsing, defaults, validation, expression
evaluation, and persistence across several APIs and solver implementations.

I propose a new configuration and experiment architecture using structured TOML,
one typed configuration model, shared staged validation, and first-class parameter
sweeps. Existing execution backends would integrate through adapters, with explicit
result/checkpoint contracts and incremental persistence migration.

**[Read the full proposal](https://github.com/ALPSim/ALPS/blob/WIP/configuration-experiments-proposal/docs/design/configuration-and-experiments.md).**

Before implementation, I would like feedback on the ownership boundaries,
compatibility commitments, sweep/replica semantics, and the first solver migration.

This issue tracks discussion and decisions; detailed specifications and revisions
live in the linked design document. The document is published on a discussion
branch without a pull request.
