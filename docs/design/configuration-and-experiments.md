# Proposal: TOML configuration, shared validation, and parameter sweeps

**Status:** Draft for discussion; no implementation or adoption is implied.

**Assessment baseline:** ALPS `master` at
[`1932d6385`](https://github.com/ALPSim/ALPS/commit/1932d63854d6e5d1a0f44c19a2b2bede45400b86).

## Overview

I propose a shared configuration pipeline for ALPS using TOML input, explicit
parameter schemas, validation before simulation startup, and typed runtime
configuration. Parameter sweeps and independent replicas would be first-class
features. Numerical results and checkpoints would continue to use HDF5.

The main change is to separate responsibilities that currently overlap:

| Responsibility | Proposed owner |
| --- | --- |
| Parse authored configuration | TOML reader or legacy input adapter |
| Describe types, defaults, and constraints | Shared framework and solver schemas |
| Expand parameter combinations and replicas | Experiment planner |
| Resolve expressions and validate scientific settings | Configuration resolver and solver validators |
| Hold effective settings for one run | Owned, typed `RunConfig` |
| Hold counters, RNG state, and measurements | Mutable solver state |
| Store numerical data and restart state | Separate result and checkpoint writers |

I recommend `toml++` for native TOML parsing and serialization. Its types should
remain internal to the reader/writer. The initial design would reuse the existing
expression engine and HDF5 archive implementation through explicit boundaries.

The first milestone should implement one complete solver path, provisionally
`spinmc`, including single runs, product and zip sweeps, validation-only planning,
resolved configuration export, and restart tests. Existing input formats and
Python APIs should remain operational during migration.

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
4. Actionable validation before expensive initialization or simulation.
5. Explicit expressions, sweep expansion, replica identity, and seed assignment.
6. A reproducible record of effective settings and resource identities.
7. Versioned result and restart contracts with a gradual compatibility path.

The initial work would not replace scientific algorithms, translate the model
and lattice description languages wholesale, replace the HDF5 wrapper, or build
a new cluster scheduler. Existing execution mechanisms should consume prepared
runs through adapters. More advanced sweep composition and conditional filtering
can follow after the basic semantics are established.

## Configuration architecture

```text
TOML file / Python mapping / legacy input adapter
                         |
         Parsed document with source information
                         |
         Apply overrides and expand experiment
                         |
         Resolve defaults and configuration expressions
                         |
         Validate and assign replica seeds/destinations
                         |
         Immutable experiment plan of typed RunConfig objects
                         |
                       Solver
                         |
          Separate result and checkpoint writers
```

The configuration core should be usable without Python, MPI initialization, or
an HDF5 archive. TOML parsing, legacy conversion, persistence, and language
bindings should be separate adapters to that core.

The parsed document should retain each value's origin: input file and position,
Python mapping path, command-line argument, or schema default. Derived values
should retain links to the inputs that produced them.

`RunConfig` should own the effective settings for one run and be immutable through
its public API. Solvers should receive typed fields rather than repeatedly
looking up and converting string keys. Dynamic model parameters can remain a
validated typed map or symbolic representation owned by the model component.
Counters, RNG state, measurements, and adaptive algorithm state belong elsewhere.

Parsing and basic validation must not open output files or construct a solver.
Scientific validation may load the selected model or lattice resource, but should
avoid simulation-sized allocations where possible. Remaining runtime invariants
must still be checked by the solver.

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
- Allowed values, bounds, and shape constraints.
- Scientific meaning and unit convention.
- Deprecated aliases and compatibility mappings.
- Whether it may change when resuming a checkpoint.

I recommend C++ schema descriptors associated with typed fields, plus named
functions for scientific constraints. Descriptors could generate reference
documentation, example files, and JSON Schema for editor assistance. The scope
should stay small rather than becoming a general-purpose reflection framework.

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

Scientific validation must be solver-specific. Examples include the relationship
between temperature and inverse temperature, update methods compatible with a
magnetic field, and admissible lattice/model combinations. A generic rule such as
requiring thermalization to be less than measurement sweeps should not be imposed
without establishing that solver's semantics.

Model/lattice parameter validation must use declarations from the selected
resource or explicit extension registrations. A fixed global allowlist would
reject legitimate user-defined models. Namespaced extension data must not become
a way to bypass validation of known solver settings.

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

Computed defaults and expressions must use final effective inputs. Canonical-name
and alias conflicts within the same layer should fail. Intentional replacement
between layers should retain the winning origin and override history.

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
cycles and undefined symbols, and resolve before solver construction when their
inputs are known. Errors should point to the authored expression. Both the
expression and effective value should be retained when applicable.

Site-, bond-, and quantum-number-dependent expressions must remain validated
symbolic objects until their evaluation context exists. Known context symbols
must be distinguished from misspelled configuration symbols.

I recommend reusing the existing expression engine through an adapter initially.
The adapter must define how namespaced references map to the existing engine.
Random expressions and disorder generation need explicit RNG ownership and
evaluation timing. Saving a configuration must never evaluate expressions again
or consume random numbers.

## Parameter sweeps and experiment planning

Here, a **parameter sweep** means executing multiple configurations. It is
distinct from `monte_carlo.sweeps`, which controls work within one simulation.

I propose supporting individual runs, explicit run lists, and parameter sweeps
with independent replicas. Planning should produce a validated persistent plan
before execution, usable by serial and existing parallel execution backends.

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
replicas = 4

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

This defines 12 parameter combinations and four replicas each: **48 runs**.

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
must specify deterministic enumeration order, although identity and seed assignment
must not depend on that order.

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
3. Resolve computed defaults and configuration expressions for each combination.
4. Validate its scientific configuration.
5. Expand replicas and assign seeds and destinations.
6. Check the complete plan for conflicts and publish it only if valid.

Expressions therefore see each combination's axis values. CLI overrides targeting
an axis must fail with an explanation. Non-axis overrides are applied before
derived values are resolved. The base input may omit fields supplied by axes;
required-field validation must account for this.

Invalid combinations must not be silently skipped. Planning should fail and
identify the coordinates and constraints involved. An eventual filtering feature
must be explicit and record excluded combinations with reasons.

### Identity and seeds

Each effective configuration should have a stable configuration identifier, and
each replica a distinct run identifier. Identity should use a versioned canonical
representation of effective scientific and sampling settings, including relevant
model/lattice resource identities. It should exclude output locations, worker
assignments, completion status, and the seed derived from that identifier.

Canonicalization must define numeric types, ordering, and expression treatment;
hashing arbitrary TOML formatting is insufficient. Full configurations should be
retained so identifiers can be checked for collisions. Execution/build metadata
should also be recorded without claiming portability of numerical trajectories.

Equivalent configurations produced more than once should be detected and rejected
initially. Additional independent sampling should be expressed through replicas.

Sweep seed assignment should derive from the experiment base seed, configuration
identifier, replica index, and a versioned derivation policy. Mapping to the RNG's
supported seed or stream representation must be specified. Duplicate effective
seed/stream assignments should be detected where applicable, and final assignments
persisted. Distinct seeds alone do not prove statistical independence.

Changing worker count, scheduling order, or backend must not change assigned
seeds. This is a seed-assignment guarantee, not a promise of identical trajectories
across solver parallelizations or platforms.

Sweep plans should reject conflicting per-run seed settings. Explicit run lists
may carry individual seeds under a separate documented policy. If automatic base
seed generation is supported, it must happen once during planning and be persisted
before execution. Legacy conversion must preserve historical resolved seeds where
required instead of silently applying the new policy.

### Inspection and persisted plans

Illustrative commands are:

```sh
alps plan experiment.toml
alps plan experiment.toml --write plan.toml
alps run --plan plan.toml
```

Planning should report axes and candidate counts, combinations, replicas, total
runs, effective settings, seeds, destinations, and validation errors. A validation
or preview operation must not create simulation output files.

The planner should calculate expansion size using overflow-checked arithmetic
before materializing combinations and enforce a documented configurable run-count
limit. Explicit value lists are sufficient initially; future range generators
must define endpoints and floating-point stepping precisely.

A persisted plan must contain resolved run definitions, resource identities,
format/schema versions, and seed policy. Execution must load that plan rather than
re-expand potentially changed authored inputs or reapply newer defaults. Plan
reading still validates integrity, compatibility, and relocated resources.

### Outputs, execution, and experiment restart

Each run should receive a unique directory:

```text
results/
  plan.toml
  runs/
    <configuration-id>/
      replica-000/
        resolved.toml
        results.h5
        checkpoint.h5
      replica-001/
        ...
```

Existing destinations require an explicit resume or replacement policy. Each
result should store experiment/configuration identifiers, replica index, sweep
coordinates, and effective settings. Analysis should read structured metadata
rather than infer it from filenames. Replicas stay distinct in storage; combining
them is an explicit statistical operation using the appropriate estimator.

Scheduling backends should consume the same prepared definitions. They must
distinguish concurrency between runs from parallelism within a run: an MPI solver
may need several ranks for one configuration.

Progress should be separate from the immutable plan, with per-run states such as
pending, running, completed, and failed. Progress publication needs coordinated
ownership; workers must not independently rewrite a shared status file. A result
should count as complete only after successful publication and validation of its
completion metadata.

Resuming an experiment should verify its plan, skip confirmed completed runs,
resume compatible incomplete checkpoints, and apply an explicit retry policy to
failed runs. Configuration, replica identity, and seed assignments must remain
unchanged. A changed experiment requires a new plan or explicit revision;
checkpoints must never be silently reassigned.

The original configuration identity remains attached to a resumed run. Permitted
restart changes, such as extending a sampling budget, should be recorded as a
validated revision, not used to derive a fresh seed or rename the existing run.

## Results, checkpoints, and provenance

| Artifact | Representation | Contract |
| --- | --- | --- |
| Authored configuration | TOML | User-maintained settings, expressions, experiment definitions |
| Resolved configuration | Generated TOML and embedded HDF5 metadata | Effective values, defaults, overrides, seeds, resource identities |
| Results | HDF5 | Typed numerical outputs and analysis metadata |
| Checkpoints | HDF5 | Complete, versioned restart state |

Serialization should store resolved values without expression evaluation or type
inference. The typed representation is authoritative; embedded text is a readable
snapshot with a defined consistency relationship, not an independently editable
second authority.

Results should record ALPS/solver and schema versions, model/lattice resource
identities or hashes, effective configuration, RNG information, and experiment
coordinates. Result and checkpoint formats need separate version identifiers:
being able to read a result does not imply being able to resume it.

Checkpoints must contain full RNG, accumulator, and solver state. Restart should
compare saved effective settings with requested changes using each field's restart
policy. New defaults must not silently reinterpret an older checkpoint. Unsupported
versions should fail with an actionable compatibility message.

Checkpoint publication should write a new generation, close it successfully, and
publish using appropriate filesystem semantics while retaining a recoverable
previous generation. Atomic visibility and durable storage are different guarantees
and must be specified for the supported environment. Loading should validate into
temporary state and leave existing simulation state unchanged on failure.

Initially preserve established HDF5 paths and encodings, including those consumed
by `/parameters` and `/simulation/results` readers. New metadata should be added
without assuming old readers understand the new schema. Independent readers should
be able to inspect supported numerical datasets without ALPS-specific decoding.

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

Python should expose experiment planning and loading of results with their sweep
coordinates and replica identities. Planning should be usable without launching
native simulations.

For MPI, the designated rank should resolve configuration and distribute the
validated representation. Failure must terminate startup coherently rather than
leave peers blocked in collectives. Execution backends must manage communicators
for parallelism within runs separately from distribution between runs. HDF5 writes
need a designated writer or explicitly supported parallel-I/O path.

## Migration and first milestone

I propose the following sequence:

1. Establish fixtures for legacy parsing, expressions, evolving parameter-list
   defaults, seeds, numeric precision, Python behavior, and checkpoints.
2. Add the configuration core and TOML reader alongside both current parameter APIs.
3. Implement the complete path for one solver, provisionally `spinmc`, using an
   adapter to existing solver interfaces where necessary.
4. Add Python preparation and legacy conversion tools with diagnostic reports.
5. Migrate remaining solvers individually, maintaining a compatibility matrix.
6. Agree on a release-based deprecation policy only after compatibility coverage.

Legacy conversion should materialize explicit resolved run definitions when the
original syntax cannot be represented faithfully as a sweep. It must report
ambiguous or unsupported transformations. TOML support must not change the meaning
of existing text files or infer the format from whether one parser happens to fail.

The first milestone should satisfy these acceptance criteria:

- A representative single-run TOML input executes `spinmc`.
- Validation-only planning creates no simulation outputs and starts no solver.
- Resolved output contains all effective settings, origins, and seeds.
- Representative legacy and TOML inputs produce equivalent effective parameters
  and scientific results under the same execution conditions.
- Product and zip sweeps expand correctly; malformed axes and conflicts fail.
- Expressions and computed defaults resolve separately per combination.
- Vector-valued parameters remain unambiguous, and numeric boundary cases are checked.
- Configuration identity and assigned seeds are stable under scheduling reordering.
- Run destinations are unique and results contain sweep/replica metadata.
- A small serial experiment can restart without rerunning completed runs or changing
  pending-run seed assignments.
- Checkpoint/resume is equivalent to uninterrupted execution within the supported
  deterministic environment; incompatible settings fail without corrupting state.
- Interrupted or failed checkpoint publication leaves a recoverable generation.
- Existing result readers and relevant Python compatibility tests continue to pass.
- Installed SDK consumers do not need TOML implementation types in solver headers.

Existing parameter, native/Python checkpoint, and HDF5 regression suites should be
reused. New tests should focus on public contracts and scientific equivalence.
Parallel execution adapters should be added with explicit multi-rank failure and
restart tests; the first milestone need not implement a new scheduling system.

## Decisions for discussion

1. Is the proposed separation of configuration, experiment planning, runtime state,
   and persistence appropriate?
2. Is `toml++` preferred, or is comment-preserving editing important enough to favor
   `toml11`?
3. Should native schema descriptors be authoritative, or should structural schemas
   be maintained externally? How should third-party solvers register extensions?
4. Are the override, expression, path, product/zip, and replica semantics suitable?
5. Which RNG seed/stream derivation policy and canonical identity representation
   should be standardized?
6. Which historical inputs, Python APIs, and checkpoint versions must remain
   supported, and for how long?
7. Is `spinmc` the best first target, and which existing execution adapter should
   follow the initial serial experiment implementation?

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

I propose a shared configuration pipeline using TOML, typed and validated run
configurations, explicit parameter-sweep planning, and versioned HDF5 result and
checkpoint contracts.

**[Read the full proposal](https://github.com/skilledwolf/ALPS/blob/WIP/configuration-experiments-proposal/docs/design/configuration-and-experiments.md).**

Before implementation, I would like feedback on the architecture, compatibility
commitments, sweep semantics, and scope of the first solver prototype.

This issue tracks discussion and decisions; detailed specifications and revisions
live in the linked design document. The document is published on a discussion
branch without a pull request.
