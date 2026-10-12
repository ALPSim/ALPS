# Roadmap: configuration and experiments

Companion to the [design proposal](configuration-and-experiments.md).
These are proposed completion gates, not release dates or one-change-per-milestone
commitments.

I propose proving one boundary at a time. The current coupling spans
[parameter parsing and persistence](../../src/alps/parameter/parameters.C),
[expression evaluation and RNG state](../../src/alps/expression/evaluator.h),
[solver initialization](../../src/apps/mc/spins/abstractspinsim.h), and
[worker construction, seeds, and restart](../../src/alps/scheduler/workertask.C).
This makes single-run correctness and restart prerequisites for experiments.

## Milestones

Complete **M0 → M1 → M2 → M3 → M4 → M5 → M6** in order. Each gate enables the
next layer. Broader coverage and input retirement have separate gates.

| Milestone | Deliver | Completion gate / next capability |
| --- | --- | --- |
| **M0 — Baseline** | Select the pilot and inventory supported files, APIs, defaults, seeds, and restart behavior. | Reproducible compatibility fixtures and a support matrix make refactoring measurable. |
| **M1 — Independent core** | Typed values, descriptors, diagnostics, immutable settings, and a resolver interface in a separate build target. | Compile/test without the monolithic ALPS library or Python/MPI/HDF5 dependencies; readers can now share the model. |
| **M2 — Input boundary** | TOML reader plus isolated historical-input importer; split shared containers/codecs from input entry points. | Supported inputs produce equivalent canonical definitions. Literal preparation works with the importer disabled; existing invocations still work in normal builds. |
| **M3 — Scientific initialization** | Extract expression resolution, graph/model setup, computed defaults, and checks from solver construction. | Preparation leaves RNG state unchanged; context-dependent checks complete before simulation. Typed settings/context can now be passed to a solver. |
| **M4 — One typed solver** | Migrate the pilot's settings and factory/task/worker boundary through the existing backend. | A serial typed run matches the baseline; each migrated field has one owner and seed transformations are accounted for. |
| **M5 — Replay and restart** | Persist prepared/effective settings, artifact roles, and complete checkpoint references; validate resume changes before loading state. | Converted plans execute, produce readable results, and resume with the importer disabled. Saving preserves RNG state; permitted updates preserve identity/seeds. |
| **M6 — Experiments and Python** | Product/zip sweeps, replicas, stable run/task mappings, owned Python snapshots, and experiment loading. | The 48-run example executes and resumes without duplicate submission or replica multiplication; analysis excludes checkpoints. This completes the first usable release. |
| **M7 — Broader coverage** | Repeat the migration for general-model solvers, additional backends, MPI, and disorder. | Each claimed capability passes equivalence and restart tests; parallel failures leave no blocked peers. Coverage grows without requiring an all-repository rewrite. |
| **M8 — Input retirement** | After conversion coverage, maintainer agreement, and a deprecation period, remove the importer from the main runtime. | Supported TOML, converted-plan, Python, result, and restart tests pass with its target excluded; migrated solver logic needs no changes. |

## Rules across milestones

- Preserve existing workflows between milestones. Unsupported features stay on
  their declared existing path; failed parsing/validation never triggers a fallback.
- Keep scientific/backend adapters separate from the removable input importer.
  Shared XML codecs and archive readers retain independent support policies.
- At M4, reconcile file-based and in-memory task setup and establish seed ownership.
  At M5, validate a whole requested settings update before applying it.
- At M6, have backends retain run IDs so an interrupted submission can be reconciled
  even if the caller failed to record the returned task mapping. Progress remains
  backend-owned.
- Use a general-model fixture during M3; `spinmc` alone cannot establish XML
  Hamiltonian coverage. General-model and parallel work may begin early, but must
  pass the M7 gates before support is claimed.
- Use targeted source splits and adapters. M8 retires an input format; removing
  remaining old parameter containers or backend formats is separate work.
