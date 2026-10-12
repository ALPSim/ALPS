# Implementation roadmap: configuration and experiments

**Status:** Proposed milestones, not implementation commitments or release dates.
Companion to the [design proposal](configuration-and-experiments.md), based on ALPS
[`master` at `1932d6385`](https://github.com/ALPSim/ALPS/commit/1932d63854d6e5d1a0f44c19a2b2bede45400b86).

I propose establishing one testable boundary at a time. Each milestone can contain
several reviewable changes; it is complete only when its exit checks pass. M6 is
the first complete user workflow. Earlier milestones must keep existing workflows
working and produce useful foundations without claiming the migration is finished.

## Dependencies that determine the order

| Current coupling | Evidence | Consequence for the plan |
| --- | --- | --- |
| Parsing, parameter storage, and persistence share classes/files | [`parameters.h`](../../src/alps/parameter/parameters.h), [`parameters.C`](../../src/alps/parameter/parameters.C), [`ngs/params.cpp`](../../src/alps/ngs/lib/params.cpp) | Split input entry points from retained containers/codecs before an importer can be optional. |
| Expression lookup uses `Parameters`; evaluator construction can seed shared disorder state | [`evaluator.h`](../../src/alps/expression/evaluator.h), [`evaluator.C`](../../src/alps/expression/evaluator.C) | Introduce explicit symbol/context/RNG ownership before promising side-effect-free preparation. |
| Graph construction, model dispatch, defaults, and checks occur during solver construction | [`montecarlo.h`](../../src/alps/scheduler/montecarlo.h), [`SpinFactory`](../../src/apps/mc/spins/spinmc_factory.C), [`AbstractSpinSim`](../../src/apps/mc/spins/abstractspinsim.h) | Extract initialization before changing solver configuration consumption. |
| Task loading, worker construction, seed updates, and resume are connected | [`task.C`](../../src/alps/scheduler/task.C), [`workertask.C`](../../src/alps/scheduler/workertask.C), [`worker.C`](../../src/alps/scheduler/worker.C) | Prove one task's execution and restart before expanding experiments. |
| Components share one `alps` library target | [`src/alps/CMakeLists.txt`](../../src/alps/CMakeLists.txt) | Start with an independent core and split only the source groups needed to enforce each boundary. |
| Python input generation assigns seeds and mutates dictionaries | [`tools.py`](../../python/pyalps/src/pyalps/tools.py) | Introduce owned snapshots through a new API while testing the existing API unchanged. |

The main dependency chain is **M0 → M1 → M2 → M3 → M4 → M5 → M6**.
M7 repeats the proven migration for further capabilities; M8 depends on coverage
of the support promise agreed in M0. General-model fixtures can start in M3 and
parallel-adapter work after M5, but neither bypasses its later acceptance checks.

## M0 — Capture the compatibility baseline

**Work:** Inventory the supported entry points: text parameter lists, XML jobs,
Python mappings, public parameter APIs, result loading, and checkpoint restart.
Record which component owns each default, conversion, seed transformation, and
state update. Select the pilot solver/backend and enumerate its promised models,
fields, and execution modes; `spinmc` with the existing scheduler is provisional.

Extend the existing [parameter tests](../../tests/parameter),
[native parameter probe](../../tests/pyalps/native_params),
[Python checkpoint contracts](../../tests/pyalps/test_checkpoint_contracts.py), and
relevant model/HDF5 tests. Include evolving defaults, empty blocks, `#clear`,
`#stop`, environment substitution, aliases, numeric boundaries, expressions,
seeds, Python mutation, and complete restart file sets. Use fixed environments,
resources, seeds, and supported topology for reproducible comparisons.

**Exit:** A checked-in support matrix and executable baseline cover the pilot's
input-to-result/restart path. Intentional behavior fixes are identified separately
from compatibility-preserving refactors. Preserve public source/link contracts
unless a separately announced change is agreed.

**Unlocks:** Every later extraction can be judged against observed behavior.

## M1 — Establish the independent configuration core

**Needs:** M0.

**Work:** Introduce a small internal CMake target for canonical typed values,
source origins, diagnostics, field descriptors, prepared settings, and immutable
effective settings. Represent deferred expressions as data with an explicit
resolution stage. Keep mutable simulation state outside these types.

Build this target without linking the monolithic `alps` library. Define a narrow
resolver interface; scientific implementations may depend on the core, while the
core does not include their `Parameters` or expression-engine headers. TOML,
Python, MPI, and archive adapters remain separate. Preserve the existing umbrella
library/package rather than restructuring every component now.

**Exit:** A standalone consumer compiles/links and exercises checked conversions,
schema composition, diagnostics, and snapshot ownership without the old library
or its HDF5/MPI/Python dependencies. CI checks the include/link boundary.

**Unlocks:** New readers and validation can develop independently of solver and
scheduler construction.

## M2 — Add TOML and isolate historical user input

**Needs:** M1 and the M0 fixtures.

**Work:** Add the pinned `toml++` reader and a separate `compat/legacy_input/`
target. Reuse the existing historical parsers, extracting their entry points from
mixed implementation files as necessary. Keep containers and XML/archive codecs
used by backends in retained shared targets. Preserve default-build public headers
and forwarding entry points, including the file constructor of `ngs::params`.

The importer expands historical parameter lists into explicit canonical run
definitions. It owns historical aliases, conversions, environment/path semantics,
and default overrides; it records changes and preserves assigned seeds. Document
normalization of accepted coercions and unused keys; dynamic scientific symbols
stay in their declared model maps. Shared schemas own current rules. Scientific
expressions remain deferred for M3; unsupported cases retain their existing execution route, selected by
declared capability rather than failed parsing/validation.

Add deterministic format dispatch and the conversion-command interface. A proposed
`ALPS_ENABLE_LEGACY_INPUT` option disables historical user-input entry points;
it does not disable backend XML/checkpoint codecs. Optional-build API limitations
must be explicit; normal builds preserve existing invocations.

**Exit:** TOML and supported historical inputs produce equivalent canonical
definitions. Literal-only preparation works with the importer disabled; there is
no link from the core to it. Existing parser fixtures and default-build downstream
API checks pass. Unsupported conversion is reported, never silently approximated.

**Unlocks:** Two input formats share one representation, with a removable import
boundary. Full scientific conversion still waits for M3–M5.

## M3 — Extract scientific resolution and initialization

**Needs:** M1–M2.

**Work:** Implement the resolver boundary using the existing expression engine's
`Evaluator` interface. Separate pure symbol lookup from RNG-dependent evaluation
and remove implicit seeding from the new preparation path. Preserve historical
seed spelling/timing through the appropriate adapter; the scheduler's
`DISORDERSEED` and Parapack's `DISORDER_SEED` are not assumed interchangeable.
Add tested namespaced-reference handling without changing historical expression
syntax globally.

Extract the pilot's graph/coupling setup, computed defaults, and scientific checks
from construction into an initialization operation returning final settings and
context. Keep dynamic model maps and site/bond/quantum symbols supported at their
declared stages. Reuse the initialized context instead of building it again inside
the solver. Add a general-model expression/substitution fixture now, since
`spinmc`'s factory does not cover general XML Hamiltonians.

**Exit:** Preparation leaves RNG state unchanged and clearly reports deferred
checks. Scientific initialization determines the same pilot defaults and rejects
the same invalid scientific combinations before simulation steps. Namespaced,
historical, cyclic, and contextual expression fixtures pass. New resolution needs
no historical input importer.

**Unlocks:** Solver constructors can consume completed settings/context instead
of interpreting an input language.

## M4 — Run the pilot through typed solver and backend boundaries

**Needs:** M3.

**Work:** Migrate `SpinFactory` model/matrix selection and the promised
`AbstractSpinSim`/`SpinSim` settings to typed consumption. Route existing inputs
through the same migrated field owners; remove duplicate defaults/checks as each
field moves. Isolate remaining graph, RNG, and scheduler `Parameters` projections
inside named adapters. Do not round-trip typed solver settings through strings.

Add a prepared-task/worker entry point to the existing backend, using its
in-memory task path where suitable. `Task`/`WorkerTask` already accept `Parameters`,
but their file and in-memory paths are not assumed behaviorally equivalent:
consolidate relevant setup/checks and adapt factory dispatch. Keep the backend's
scheduling decisions. Establish exactly one owner for task-to-worker seed handling;
account for the existing `WorkerTask` seed increments before any replica feature.

**Exit:** A serial single run reaches the actual simulation loop using typed
settings, on the existing backend, and agrees with the M0 fixture. The promised
field set has one owner; remaining legacy consumers are explicitly listed. No
TOML parser or user-input importer is needed by the migrated solver. Existing
unsupported models/modes continue through their original entry points.

**Unlocks:** Persistence can describe settings actually used by a running solver.
This is an internal execution milestone; complete supported restart is M5.

## M5 — Make single-run persistence and restart independent of input

**Needs:** M4.

**Work:** Persist canonical prepared definitions and final effective settings with
run IDs, destinations, versions, resources, origins, seeds, and artifact roles.
Historical conversion must leave no executable dependency on the importer or opaque legacy input blob;
required scientific resolvers are recorded separately.

Give the pilot backend metadata inspection and restart-validation entry points.
Check requested changes against saved effective settings before loading worker
state; do not recompute defaults or reseed while interpreting a checkpoint.
Inventory all files required for restart. Retain supported backend state formats
and document their publication/recovery limits.

Adapt the migrated `Worker::set_parameters()` path to validate a complete requested
update before applying it. Accepted budget changes produce a new immutable settings
snapshot and recorded update; they do not reset counters or RNG state. Historical
callbacks remain confined to unmigrated paths.

Keep `/simulation/results` and checked analysis projections such as `/parameters/T`.
Their writers belong to persistence adapters, independently of historical input.
For the migrated path, replace parameter serialization that evaluates expressions
or infers types; verify that saving does not advance RNG state. Preserve historical
reader support separately.

**Exit:** Prepare, execute, analyze, and resume the pilot with the importer disabled,
including a plan converted while it was enabled. Editing the original input cannot
change replay. Incompatible requests fail before state loading; permitted budget
extensions retain identity/seeds. Fixed-topology restart preserves the supported
scientific behavior. Archive readers pass their existing compatibility fixtures.

**Unlocks:** An experiment can safely refer to independently replayable runs.

## M6 — Add sweeps, experiment execution, and Python

**Needs:** M5.

**Work:** Extend single runs to explicit lists, product/zip expansion, count limits,
per-point defaults, and sampling replicas. Persist run definitions, coordinates,
task seeds, and destinations before submission. Have the backend retain run IDs
when creating tasks, then record returned task mappings. This permits reconciliation
if submission succeeds but recording the mapping is interrupted. Reconcile repeated
submission/resume against backend state; progress remains backend-owned. Define
and test the recovery behavior around each submission interruption point.

Fix the mapping between each planned sampling replica and backend workers/chains;
prevent backend expansion from multiplying replicas. Begin with the validated
serial pilot topology, explicitly rejecting unsupported disorder/parallel layouts.
New plan/run commands and Python mapping/file APIs share native preparation and
owned snapshots. Add artifact-aware experiment loading without checkpoint files
appearing as results; keep the old mutable Python API working.

**Exit:** The proposal's 12-point/four-replica example executes 48 intended runs.
Zip/error/limit cases pass; scheduling order and permitted resume preserve IDs and
task seeds. Interruption/retry does not duplicate submission under the documented
recovery contract. Python analysis preserves coordinates and replica identities.
Run both compatibility-enabled and importer-disabled end-to-end CI.

**Unlocks:** The first usable release, with typed solver consumption, sweeps,
analysis, conversion, and restart. Publish exact solver/mode coverage with it.

## M7 — Extend coverage through the established boundaries

**Needs:** M6 for release; general-model fixtures and parallel-adapter development
may start earlier as described above.

**Work:** Apply the same sequence to further solver/model combinations, including
a general XML Hamiltonian consumer, and then additional backend capabilities.
For Parapack, explicitly map task/clone/worker identities, separate disorder and
sampling seeds/aggregation, and account for HDF5 plus XDR checkpoint sets. For MPI,
test settings distribution and coordinated initialization failure as well as
simulation/restart. Preserve backend-specific stream behavior rather than promising
invariance across process counts.

**Exit:** Each newly claimed capability has old/new equivalence, initialization,
seed, result, and restart fixtures. Parallel failures leave no blocked peers.
Remove internal `Parameters` dependencies only as their consumers migrate; require
a dependency audit and downstream tests for each removal. Update the support matrix
after each solver/backend, not only after an all-repository rewrite.

**Unlocks:** Evidence for expanding default support and considering input retirement.
Individual checkpoint backends can separately adopt complete generation publication
and transactional loading after their state contracts are covered.

## M8 — Retire historical input without retiring archives

**Needs:** Conversion/execution coverage for the M0 support promise, explicit
maintainer agreement, and an announced deprecation period. M6 alone is insufficient.

**Work:** Finish conversion diagnostics and user migration instructions. Default
builds retain old inputs throughout the support window. After it ends, remove the
input registration, compatibility target, and forwarding APIs according to the
announced source/binary compatibility policy. The importer may remain available in
an optional conversion tool; it must not be a dependency of the main runtime.

**Exit:** A clean build with the importer target physically excluded passes the
supported TOML, converted-plan, Python, result, and restart suites. Migrated solver
logic requires no changes. Remaining shared XML/expression/backend codecs and
historical archive readers retain their independent support policies.

**Unlocks:** Subsequent removal of individual old parameter containers or backend
formats can proceed on their own evidence and schedule.

## Delivery rules

- Give each extraction a behavior-preserving change and focused checks before
  routing a public entry point through it. Keep all promised legacy routes working
  between milestones; never fall back to another parser after an input error.
- Avoid a plugin framework, a versioned interpreter for every historical release,
  or a repository-wide directory move. Use concrete adapters and targeted source
  splits with explicit dependency direction.
- M1–M6 may merge incrementally, but the user-facing release claim requires all M6
  gates and the [proposal's acceptance checklist](configuration-and-experiments.md#migration-and-acceptance).
  Version support, scientific equivalence tolerances, and recovery guarantees must
  be stated; unexplained differences block expansion of the supported set.
- Follow the [contribution and provenance policy](../../CONTRIBUTING.md#provenance-and-scientific-credit).
  Moving existing code preserves its notices; these documents add no implementation.
