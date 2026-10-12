# Proposal: TOML configuration, shared validation, and parameter sweeps

**Status:** Draft for [discussion in issue #214](https://github.com/ALPSim/ALPS/issues/214).
Based on [`master` at `1932d6385`](https://github.com/ALPSim/ALPS/commit/1932d63854d6e5d1a0f44c19a2b2bede45400b86).

I propose structured TOML input, one typed configuration model, shared validation,
and first-class parameter sweeps. Existing files keep working during migration.
Existing backends retain scheduling and restart responsibilities; HDF5 remains the
result format. The [roadmap](configuration-migration-roadmap.md) breaks the work
into independently testable milestones.

Today, parameter interpretation is spread across readers, solver constructors,
Python helpers, and persistence. The change gives each decision one owner and makes
it possible to inspect what will run before starting a simulation.

## Input and ownership

Use [toml++](https://github.com/marzer/tomlplusplus) with a versioned TOML 1.0
profile. Keep parser types outside public solver APIs; pin the dependency and
support offline builds. This illustrative input defines 12 parameter points with
four sampling replicas each: **48 runs**.

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
Names and interfaces remain open for discussion.

| Component | Responsibility |
| --- | --- |
| Configuration core and solver schemas | Typed values, defaults, constraints, and diagnostics |
| Historical-input importer | Old syntax and behavior translated into canonical definitions |
| Experiment planner | Parameter points, replicas, run IDs, seeds, and destinations |
| Scientific initialization | Context-dependent defaults/checks and final settings |
| Existing execution backend | Scheduling, progress, retries, and continuation |
| Persistence adapters | Results, checkpoints, and compatibility with existing readers |

The core is independent of Python, MPI, HDF5, and old parameter containers.
Prepared settings are immutable; explicitly deferred values are resolved during
scientific initialization. Final settings stay separate from mutable solver state.
Each migrated field has one owner for defaults and checks.

## Validation and expressions

Use C++ field descriptors and scientific validators shared by TOML and Python.
Validation happens at three stages:

- **Preparation:** types, required fields, available constraints, and sweep errors.
- **Initialization:** checks/defaults requiring a graph, model, or algorithm context.
- **Runtime:** invariants involving evolving simulation state.

Defaults apply only to absent fields. Conversions must detect overflow, truncation,
and precision loss. Reject unknown fields where the schema is complete; preserve
dynamic scientific symbols where model/lattice definitions require them. Errors
identify the field, source, and affected run, including checks still deferred.

Expressions are explicit, for example `J = { expr = "2 * constants.J0" }`.
Adapt the existing evaluator for namespaced references and contextual symbols;
detect undefined references and cycles. Preparation and serialization must not
change random-number state. Random expressions require explicit evaluation timing
and RNG ownership.

Record value origins and overrides. Resolve file paths relative to their source
file. New Python APIs use owned snapshots, so subsequent dictionary/array mutations
cannot change prepared runs. Existing mutable APIs retain their behavior during
migration.

## Parameter sweeps

Support single runs, explicit run lists, and two sweep modes: **product** takes all
combinations; **zip** pairs entries from equally sized axes. Ordinary arrays remain
parameter values. A parameter sweep is distinct from `monte_carlo.sweeps`, the
work performed within one run.

Expand in a defined order, resolve defaults/expressions per point, then validate.
Check counts before allocating runs. Reject empty axes, conflicting overrides,
invalid candidates, and unsupported targets; never silently discard invalid points.
Initially restrict axes to scientific/solver settings and keep explicit lists and
sweeps separate.

Preview shows run counts, settings, seeds, destinations, and deferred checks without
simulation outputs. Save the prepared experiment before execution so edited inputs
or changed defaults cannot silently alter replay.

Keep parameter points, disorder realizations, sampling replicas, and workers
separate. Sampling replicas do not imply new disorder realizations. Persist simple
IDs and task seeds once; preserve them across scheduling and supported restart.
Backend-specific worker streams may depend on topology. Adapters must prevent
backend-created chains from multiplying planned replicas and reject unsupported
layouts before submission.

The backend remains the sole authority for progress. Map run IDs to its tasks and
reconcile interrupted submission/resume to avoid duplicate execution.

## Results and restart

Keep three distinct records: authored input, prepared runs, and effective settings
actually used. Results record run identities, coordinates, replicas, versions,
resources, and seeds. Saving preserves types and never evaluates expressions.

Artifacts have explicit roles. The new loader discovers declared results without
mistaking checkpoints for measurements. Preserve compatible numerical datasets and
analysis metadata, including checked projections such as `/parameters/T`. Output
destinations must be unique; reusing them requires explicit resume or replacement.

A checkpoint must identify the complete restart state, including RNGs, measurements,
and solver state. Existing backends may require multiple HDF5/XDR files. Validate
requested changes against saved effective settings before loading state; supported
budget extensions retain identity and seeds. Result readability and checkpoint
resumability have separate versioned contracts.

Test the selected backend's existing restart behavior and document its recovery
limits. More robust checkpoint publication and uniform storage can follow per backend.

## Existing inputs and eventual removal

Keep old input support in a separate `compat/legacy_input/` build target. Reuse its
parsers and preserve historical defaults, expressions, environment/path behavior,
run ordering, and seed assignments. Convert parameter lists to explicit runs when
that preserves their meaning. Unsupported cases retain their existing execution
path until migrated; never retry another format after a parse error.

The importer depends on the canonical model; the core and migrated solvers never
depend on the importer. Historical defaults become explicit values where possible.
Deferred scientific resolution belongs to a supported resolver. Converted plans
must execute with the importer disabled.

Shared expression/model/lattice support, backend serialization, and historical
archive readers have separate lifetimes. Removing old input support must not remove
them. Existing commands and public APIs remain supported through the announced
compatibility period; conversion tooling reuses the importer and leaves originals
intact. Retirement requires demonstrated coverage and maintainer agreement.

## First release and open decisions

The first usable release, M6 in the [roadmap](configuration-migration-roadmap.md),
must demonstrate a solver consuming typed settings, sweeps, Python analysis, and
restart on an existing backend. Provisionally use clean-model `spinmc` with serial
execution. Verify historical/TOML equivalence, existing reader/Python contracts,
and execution of converted plans with the importer disabled. Publish exact coverage.

Broader solvers, MPI, and disorder support follow with their own tests. Advanced
sweep modes, generated editor schemas, content-derived IDs, full provenance graphs,
and cross-backend checkpoint migration are deferred.

Decisions needed: agree on the ownership boundaries, TOML/schema approach, first
solver/backend, and compatibility scope/support period. Detailed API and storage
choices belong to the implementation milestones.

Implementation follows the [contribution and provenance policy](../../CONTRIBUTING.md#provenance-and-scientific-credit).
