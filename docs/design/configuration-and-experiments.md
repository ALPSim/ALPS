# Proposal: TOML configuration, shared validation, and parameter sweeps

**Status:** Draft for discussion; no implementation or adoption is implied.
**Baseline:** ALPS `master` at
[`1932d6385`](https://github.com/ALPSim/ALPS/commit/1932d63854d6e5d1a0f44c19a2b2bede45400b86).

## Overview

I propose a new configuration and experiment architecture for ALPS: structured
TOML, one typed configuration model, shared staged validation, and first-class
parameter sweeps. Existing solvers and execution backends would migrate through
adapters. HDF5 remains the result format and the target for migrated checkpoints.

**The experiment model defines what should run; the execution backend owns its
progress.** Existing parameter containers become compatibility boundaries. The
first usable release must migrate one solver to typed configuration, provisionally
`spinmc`, and demonstrate sweeps, analysis, and restart through its existing backend.

I recommend `toml++`, retaining the current expression engine and HDF5 wrapper
initially. Content-derived IDs, uniform checkpoint storage, and advanced plan
revision tools can follow later. All new APIs and examples below are illustrative.

## Contents

- [Why change](#why-change)
- [Architecture and ownership](#architecture-and-ownership)
- [Configuration contract](#configuration-contract)
- [Sweeps and execution](#sweeps-and-execution)
- [Results and restart](#results-and-restart)
- [Migration and acceptance](#migration-and-acceptance)
- [Open decisions](#open-decisions)
- [References and provenance](#references-and-provenance)
- [Short issue text](#short-issue-text)

## Why change

Supported parameters, defaults, and effective settings are difficult to discover
because their interpretation is distributed across ALPS:

- [`Parameters`](../../src/alps/parameter/parameters.h) stores text;
  [`params`](../../src/alps/ngs/params.hpp) adds typed values but still
  [uses the legacy text parser](../../src/alps/ngs/lib/params.cpp).
- Solver constructors such as
  [`AbstractSpinSim`](../../src/apps/mc/spins/abstractspinsim.h) own defaults and
  scientific checks, while [`Parameters::save()`](../../src/alps/parameter/parameters.C)
  evaluates expressions and chooses numeric representations during serialization.
- [`writeInputFiles()`](../../python/pyalps/src/pyalps/tools.py) combines seed
  assignment, dictionary mutation, job generation, and XML writing.
- Existing contracts must survive migration: evolving defaults in
  [`ParameterListParser`](../../src/alps/parameter/parameterlist_p.h),
  [Python mutation and checkpoint behavior](../../tests/pyalps/test_checkpoint_contracts.py),
  and established [result layouts](../../python/pyalps/src/pyalps/load.py).

Execution and persistence already have owners. [Parapack tasks](../../src/alps/parapack/job.C)
manage clones, continuation, and aggregation; its
[seed initialization](../../src/alps/parapack/clone_info.C) separates worker and
disorder seeds. [Clone persistence](../../src/alps/parapack/clone.C) can combine
HDF5 metadata with XDR worker state. The new architecture must integrate these
capabilities explicitly.

## Architecture and ownership

| Component | Owns |
| --- | --- |
| Configuration core and solver schemas | Typed values, defaults, conversions, constraints, and source locations |
| Optional historical-input importer | Old syntax, historical names/conversions/defaults, and conversion to canonical run definitions |
| Experiment planner | Parameter points, run IDs, coordinates, replica definitions, task seeds, and output destinations |
| Scientific initialization | Context-dependent defaults and checks; final immutable `RunConfig` |
| Existing execution backend | Scheduling, progress, completion, retries, and continuation |
| Solver state | Mutable counters, RNG state, measurements, and adaptive state |
| Result/checkpoint adapters | Artifact roles, versioned persistence contracts, and compatibility |

```text
TOML / Python mapping / legacy adapter
  -> typed document with source locations
  -> overrides, sweep expansion, preparation checks
  -> persistent definition of prepared runs
  -> execution backend -> scientific initialization
  -> solver using RunConfig -> results/checkpoints
```

The core should require neither Python, MPI initialization, nor an HDF5 archive.
Prepared settings are immutable; context-dependent values stay explicitly deferred
until initialization returns a final `RunConfig` and scientific context. Effective
output records those decisions without rewriting authored intent. Dynamic model
parameters may remain typed maps or symbolic objects.
Value origins identify a file position, Python mapping path, CLI argument, or
schema default.

Each migrated field has one owner for defaults and checks. Adapters may translate
settings for remaining legacy dependencies, but must have a migration path: the new
model should replace the existing parameter APIs progressively. A TOML-to-legacy
converter alone does not complete the first usable release.

Input compatibility lives in a separate target under `compat/legacy_input/`, which
depends on the canonical model; the core never depends on the importer. Historical
defaults become explicit imported values where possible. Deferred scientific
resolution belongs to a supported resolver, so converted plans can execute without
the importer. No general `legacy_mode` flag crosses into solver logic.

Shared expression/model/lattice code, backend serialization, and historical archive
readers have separate lifetimes. Existing `Parameters` consumers remain behind
adapters until migrated. Disabling the user-input importer does not remove shared
XML codecs or backend checkpoint readers. The [implementation roadmap](configuration-migration-roadmap.md)
specifies the dependency boundaries and their removal tests.

Python uses the same native preparation rules and receives owned snapshots; later
dictionary/array mutations cannot change prepared runs. Existing mutable
`ngs.params` behavior remains in compatibility APIs. Arrays, nested metadata, and
exact integers need checked mappings; unsupported values must not be silently
stringified or narrowed. User annotations remain separate from solver settings.

## Configuration contract

### Format and dependencies

Use a documented TOML 1.0 profile, with experimental extensions disabled and unknown
format versions rejected. [`toml++`](https://github.com/marzer/tomlplusplus) fits
C++17 and provides source locations; keep its types out of public solver headers.
Pin the dependency, support offline/system-package builds, and retain its notices.
[`toml11`](https://github.com/ToruNiina/toml11) is an alternative if editing files
while preserving comments becomes essential. Initially, write separate snapshots
and leave authored files intact.

This complete example defines 12 parameter points and four sampling replicas each:
**48 runs**.

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
temperature = 0.5
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

For a single run, omit `experiment` and `sweep`, and set `monte_carlo.seed`.
Individual output paths may also be specified. Existing XML model/lattice libraries
remain accessible through references.

Keep scientific symbols such as `J` and `L`. Document legacy mappings such as
`monte_carlo.temperature` -> `T` and `monte_carlo.sweeps` -> `SWEEPS`, including
analysis metadata. Historical input aliases belong to the importer; backend and
output projections must work independently of it. Reject ambiguous mappings; do
not flatten away namespaces.

The value profile must specify:

- Signed 64-bit TOML integers. Larger historical/Python values require explicit
  compatibility representation or checked rejection, never a floating-point detour.
- Dates/times only in fields that admit them; complex values through an explicit
  representation such as `{ real = 1.0, imag = 0.5 }`.
- Typed dataset references for large arrays rather than embedding them in TOML.

### Schemas and validation

I recommend C++ field descriptors plus named scientific validators. A descriptor
declares the path/type, whether required, literal or computed default, bounds/shapes,
meaning/units, owning component, resolution stage, and restart policy.
Common and solver-specific schemas compose; generated documentation/editor schemas
can follow. External structural schemas remain an alternative, but would need a
specified dialect and default policy plus the same scientific validators.

| Stage | Checks and responsibilities |
| --- | --- |
| Preparation | Syntax, known fields, required values, checked conversions, and available cross-field constraints; no solver construction or simulation outputs |
| Initialization | Defaults/checks requiring graph, model, or algorithm context; complete before simulation steps |
| Runtime | Invariants involving evolving simulation state |

For example, [`spinmc`'s automatic update choice](../../src/apps/mc/spins/abstractspinsim.h)
depends on graph frustration and couplings. Share that rule between initialization
and optional thorough validation instead of duplicating it in preflight checks.
Temperature/beta relationships and field/update compatibility also belong to their
solver; do not impose universal rules such as thermalization < measurement sweeps.

Defaults apply only to absent fields. Distinguish booleans, integers, fractions,
overflow, and precision loss: `sweeps = 1000.5` must fail. Any integer-to-real
conversion policy must be explicit; infinity is allowed only where meaningful.

Reject unknown fields where the owning component has a complete schema. XML
defaults are not complete model schemas: expressions and substitutions such as
`J#` can introduce symbols. Legacy model adapters must declare validation coverage,
preserve symbol resolution, and diagnose unresolved symbols at the proper stage.
Extensions must not bypass checks on known controls.

Collect independent errors with source location, field path, expected/actual value,
stable error code, and run coordinates. Suppress cascading errors and identify
checks still deferred. For example:

```text
run.toml:18:10: monte_carlo.sweeps
  expected a positive integer; received 1000.5
run.toml:19:1: unknown field "thermalisation"
  did you mean "thermalization"?
```

### Overrides, paths, and expressions

For single runs and explicit run lists, precedence is:

```text
schema defaults -> base input -> per-run overrides -> explicit CLI overrides
```

Reject same-layer alias/canonical-name conflicts; record winning origins and
applied overrides. Derived values use final inputs at their declared stage.
File paths resolve relative to their source file, CLI paths to the invocation
directory. Relocation must be explicit and check resource identity. New inputs
have no implicit environment substitution; any future explicit substitution must
record effective values. Legacy adapters preserve existing behavior.

Expressions are explicit rather than inferred from strings:

```toml
[constants]
J0 = 1.0
[model.parameters]
J = { expr = "2 * constants.J0" }
```

Define namespaces, detect cycles/undefined references, preserve expressions and
resolved values, and report failures at the authored location. Site/bond/quantum
symbols remain symbolic until their context exists; distinguish allowed context
symbols from unknown references. Reuse the existing evaluator
through an adapter: its [identifier grammar](../../src/alps/parser/parser.h) does
not accept dotted paths, so namespaced references need tested resolver/parser
integration. Random expressions and disorder require explicit RNG ownership and
evaluation timing. Serialization must never evaluate expressions or consume RNG state.

## Sweeps and execution

A **parameter sweep** runs different configurations; `monte_carlo.sweeps` controls
work within one run. Support single runs, explicit run lists, and the following
expansion modes:

| Mode | Behavior |
| --- | --- |
| `product` | Cartesian product of axis values |
| `zip` | Corresponding entries, with equal-length axes and no broadcasting |

For example, zipping `L = [32, 64, 128]` with
`sweeps = [100_000, 200_000, 400_000]` produces three points. Ordinary arrays remain
parameter values; each candidate for a vector-valued axis is a complete vector.

Axis targets are declared scientific/solver fields, addressed by configuration
paths with quoted segments where needed. Reject empty axes, duplicate/unknown
targets, and invalid candidate types/shapes. An axis may supply a required field
missing from the base. Initially exclude solver selection, versions, replica counts,
seeds, and output destinations as targets; reject combining a sweep with an explicit
run list. Grouped modes, filtering, and ranges are later extensions; ranges would
need precise endpoint and floating-point stepping rules.

### Preparation and inspection

1. Apply base and experiment-wide overrides. Reject CLI overrides of axis targets.
2. Expand axes in a specified deterministic order, replacing their base values.
3. Resolve preparation-time defaults/expressions per point, then validate inputs.
4. Expand sampling replicas; assign IDs, task seeds, and destinations.
5. Check conflicts/limits and persist the definition, including deferred checks.

Calculate run counts with overflow checks and a configurable limit before
materializing combinations. Invalid points must not disappear: preparation fails
with coordinates and reasons; later initialization failures become backend-reported
failed runs. Future filtering must explicitly record exclusions.

Proposed commands:

```sh
alps plan experiment.toml
alps plan experiment.toml --write plan.toml
alps run --plan plan.toml
```

Preview reports axes/counts, settings, seeds, destinations, errors, and deferred
values/checks without creating simulation outputs. Optional thorough validation
initializes scientific context without simulation steps.

The persisted definition contains prepared settings, deferred values, IDs,
resource identities, schema/resolver versions, and assigned seeds. Execute it
without re-expanding edited inputs. Changed default semantics require compatible
resolution, explicit migration, or rejection. Final initialization decisions are
saved as the effective configuration used for restart.

### Identity, replicas, and seeds

| Identity | Meaning |
| --- | --- |
| Parameter point | One combination of scientific and sampling settings |
| Disorder realization | A particular random physical instance |
| Sampling replica | A separate chain for that physical instance |
| Worker/rank | An execution participant within a run |

Persist an experiment ID, point IDs, distinct run IDs, coordinates, and replica
indices. Explicit or simply generated IDs suffice; configuration hashes may later
provide provenance but do not define identity/seeds. Warn on simple duplicate
points, allowing deliberate repeats with distinct IDs; no symbolic-equivalence
analysis is needed.

Assign task seeds once from the recorded base seed and preparation order using a
documented policy; explicit run lists may supply them directly. Check backend seed
ranges and conflicting declarations. Automatic base-seed generation happens once
and is saved. Scheduling/restart preserves persisted IDs and seeds; re-preparing
edited/reordered inputs need not. Legacy imports retain historical assignments and
backend seed policies.

Sampling replicas must not imply new disorder realizations. Disorder-capable
adapters must specify realization seeds and sharing across points. Record sampling
and disorder contributions separately when aggregating, reusing existing clone
evaluators where appropriate. The initial clean-model `spinmc` path supports
sampling replicas and explicitly rejects unsupported disorder requests.

Adapters document task/clone/worker mappings and prevent planned replicas from
multiplying backend-created chains accidentally. Reject unsupported layouts before
submission. The guarantee is stable task-level seed assignment; existing worker
streams may depend on topology. Neither stream invariance across backends/process
counts nor statistical independence from distinct seeds is assumed.

### Backend integration

Map persisted run IDs to backend task/clone IDs, reconciling submission and resume
to avoid duplicate execution. Expose backend progress through a common view without
creating another writable status database. The backend owns completion, retries,
and continuation; document its recovery limits.

Backends also own communicators and scheduling between/within runs. Under MPI,
prepare/distribute typed settings on the designated rank and coordinate
initialization checks across participants. Fail coherently before simulation steps
so peers cannot remain blocked. HDF5 writes require a designated writer or an
explicitly supported parallel-I/O path.

## Results and restart

| Artifact | Representation and authority |
| --- | --- |
| Authored input | TOML; user-maintained settings and expressions |
| Prepared experiment | Persistent typed definition with TOML export; runs, assignments, deferred checks |
| Effective configuration | TOML snapshot and HDF5 metadata; final settings actually used |
| Results | HDF5; numerical outputs, analysis metadata, and contributing run identities |
| Checkpoints | Backend-specific during migration, HDF5 for migrated implementations; complete restart state |

Preserve typed values/deferred expressions without type inference on save. Record
the relationship between prepared intent and effective settings; readable exports
cannot become conflicting restart authorities. Keep source locations, defaults,
overrides, and resource identities. Full provenance dependency graphs and semantic
hashing are optional later.

Results record experiment/point/run IDs, coordinates, sampling/disorder identities,
effective configuration, ALPS/solver/schema versions, resources, and RNG information.
Artifacts have explicit roles and backend/format identifiers; result readability
and checkpoint resumability are separate versioned contracts.

Output destinations must be unique and require explicit resume/replacement when
already occupied. Layout is optional: use per-run directories or backend filenames,
recording artifact locations and roles. The new loader enumerates declared results
and checks roles. Legacy [`getResultFiles()`](../../python/pyalps/src/pyalps/tools.py)
may match all HDF5 files, including checkpoints, so its heuristic stays confined to
legacy workflows. Preserve numerical datasets such as `/simulation/results` and
generate checked compatibility projections such as `/parameters/T`; existing
per-file readers and independent HDF5 readers can then remain useful.

A checkpoint includes all RNG, accumulator, and solver state, potentially across
HDF5/XDR files. Adapters identify the complete set. Before loading state, validate
requested changes against saved effective settings, shared field policies, and
backend capabilities. Unsupported versions or changes fail clearly; new defaults
must not reinterpret saved state. Allowed budget extensions retain IDs/seeds and
record a configuration update.

Initially preserve and test the selected backend's restart behavior, documenting
publication/recovery limits. Migrated implementations should publish closed,
complete generations, retain a recoverable predecessor, and validate temporary
state before replacing live state. Specify atomic visibility and durability
separately. These are incremental persistence milestones; existing multi-file
writes must not be described as already transactional.

## Migration and acceptance

I propose dependency-gated milestones, detailed in the
[implementation roadmap](configuration-migration-roadmap.md). Each has source-level
scope, completion checks, and a specific next capability it enables.

| Milestone | Unlocks |
| --- | --- |
| M0: Capture existing behavior and dependencies | Refactoring against a reproducible compatibility baseline |
| M1: Introduce the independent typed core | Configuration work outside the monolithic library |
| M2: Add TOML and isolate historical input | Both formats produce canonical definitions |
| M3: Separate scientific resolution from construction | Shared validation and initialization without simulation steps |
| M4: Migrate one solver and its task/worker boundary | A typed single run on the existing backend |
| M5: Persist definitions, results, and restart decisions | Replay/resume without the input importer |
| M6: Add experiments and Python integration | The first usable release: sweeps, analysis, and restart |
| M7: Extend solver/backend coverage | Evidence for broader compatibility and parallel/disorder support |
| M8: Retire historical input separately | Remove its target without changing migrated solver logic |

M0–M5 are independently reviewable foundations; M6 is the first complete user
workflow. General-model coverage and retirement require their own later gates.

Legacy conversion must preserve evolving parameter-list defaults and materialize
explicit runs when no faithful sweep exists. Report ambiguous/unsupported cases;
retain their existing execution path until supported. Keep existing invocations
working; do not choose formats by parser-failure fallback. A conversion command
reuses the importer and reports changes without rewriting the original file.

The first usable release (M6) is accepted when:

| Area | Required evidence |
| --- | --- |
| Core | One solver consumes typed configuration; migrated defaults/checks have one owner; public solver headers hide TOML types |
| Validation | Input preview has no simulation side effects and reports deferred checks; initialization gates execution; numeric boundaries, vectors, namespaced/legacy expressions, and per-point defaults are checked |
| Sweeps | Product/zip counts and errors are correct; IDs/task seeds persist across execution order and permitted resume changes; repeats have distinct IDs; replica/backend mapping and unsupported disorder/topologies are explicit |
| Persistence | Prepared and effective settings are recorded; destinations are unique; results have role/coordinate/replica metadata; new discovery excludes checkpoints |
| Execution | The chosen backend executes/resumes an experiment without duplicate submission or a second progress database; a deterministic restart fixture preserves supported scientific behavior and rejects incompatible settings before state loading |
| Compatibility | Historical invocations and representative legacy/TOML results agree; existing per-file readers and Python contracts pass; converted plans execute with the input importer disabled; recovery limits and supported features are documented |

Reuse existing parameter, native/Python checkpoint, and HDF5 suites. Add multi-rank
failure/restart tests as parallel adapters migrate. Python's file and mapping entry
points must share native rules; file paths needing identical diagnostics use the
native reader. Optional utilities use `tomllib` on Python 3.11+ or Tomli on older
supported versions, Tomli-W for generation, and TOML Kit for formatting-preserving
editing. Alternative parser paths must pass the same format-profile tests.

**Deferred:** generated editor schemas, full provenance graphs, content-derived
identity, grouped/range/filtered sweeps, advanced plan revisions, cross-backend
checkpoint migration, generation retention, and transactional loading across all
backends. Replacing scientific algorithms, the scheduler, model/lattice description
languages wholesale, or the HDF5 wrapper is outside the initial implementation.

## Open decisions

1. Is the ownership split between configuration/experiments and execution suitable?
2. Is `toml++` preferred, or does comment-preserving editing favor `toml11`?
3. Should C++ descriptors be authoritative, and how should solver/model extensions
   and staged validators register?
4. Are the names, legacy mappings, expressions, paths, and sweep semantics suitable?
5. Are point/disorder/replica/worker identities sufficient, and which backend
   task-seed policy should the first adapter use?
6. Which historical inputs, Python APIs, and checkpoints need support, for how long?
7. Is `spinmc` the right first target, and which general-model solver/backend follows?

## References and provenance

This proposal is based on the linked ALPS source at the stated baseline. Examples
describe proposed interfaces; this branch adds no dependency or copied external
implementation.

- [TOML 1.0](https://toml.io/en/v1.0.0),
  [toml++ source locations](https://marzer.github.io/tomlplusplus/structtoml_1_1source__region.html),
  [Python TOML utilities](https://docs.python.org/3.11/library/tomllib.html).
- [JSON Schema validator](https://github.com/pboettch/json-schema-validator), an
  option if external structural schemas are selected; its dialect/default behavior
  must be specified.
- [HighFive](https://github.com/highfive-devs/highfive), a separate possible archive
  modernization with its own compatibility tests.

Implementation must follow the [contribution policy](../../CONTRIBUTING.md#provenance-and-scientific-credit):
preserve provenance/license notices and maintain relevant scientific credit in
`CITATION.cff` and `CITATIONS.yaml`. This document changes neither scientific
implementations nor citation ownership.

## Short issue text

**Title:** Proposal: TOML configuration, shared validation, and parameter sweeps

ALPS distributes parameter parsing, defaults, validation, expressions, and
persistence across several APIs and solvers. I propose structured TOML, one typed
configuration model, shared staged validation, and first-class sweeps, with
incremental migration through existing execution backends.

**[Read the full proposal](https://github.com/ALPSim/ALPS/blob/WIP/configuration-experiments-proposal/docs/design/configuration-and-experiments.md).**

I would like feedback on ownership, compatibility, sweep/replica semantics, and
the first solver migration. Discussion and decisions belong in this issue; the
detailed design lives on the linked branch, without a pull request.
