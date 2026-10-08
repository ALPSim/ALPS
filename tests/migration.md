# Native compatibility references

The GoogleTest suites preserve historical scenarios while replacing printed
numeric output with direct assertions. This file records reference provenance
and compatibility decisions needed when changing implementations. Run commands
and test conventions are in [the testing guide](README.md).

## Reference data and tolerances

| Area | Reference contract |
| --- | --- |
| Matrix algorithms | Preserve the original scalar types, mixed-type pairs and shape combinations. The former Boost percentage tolerance of `1e-6` is a relative bound of `1e-8`; identity residuals use an absolute `1e-6` bound. Seed 3 makes inputs independent of execution order. |
| Matrix archives | Retain the immutable `matrix_deprecated_hdf5_format_test` input, all stored entries, dimensions and capacity. Decimal comparisons allow the former fixture's six-significant-digit rounding. Current-writer tests separately check disk axes, complex layout and empty matrices. |
| Containers | Compare modifier operations and seeded operation traces against standard containers. Non-POD lifetime checks track live addresses without reading uninitialized objects. |
| Alea | Numerical tables come from the pre-migration stdout fixtures, not newly generated implementation output. Tolerances reflect historical print precision: six significant digits in the observable tests and ten in `mcdata2`. They establish deterministic compatibility, not statistical acceptance intervals. |
| Parameters and expressions | Preserve grammar, ordering, substitution, seeded evaluation and canonical symbolic strings. Analytic values supplement decimal references. HDF5 and XDR cases check decoded values as well as serialization. |
| Parapack | Preserve schedules, worker metadata, merging, signed logarithmic arithmetic and seeded Wang–Landau results. Centered regression calculations provide independent references. The merged mean's absolute `1e-8` bound covers roundoff. |
| XML | Keep exact XML and symbolic-expression fixtures where spelling, ordering or schema is the compatibility contract. `StreamFixture` restores streams during exceptions and compares output after serialization objects are destroyed. |

Do not regenerate historical fixtures from a replacement implementation merely
to make a test pass. The frozen Python checkpoint has its serializer revision
and reproduction instructions beside it in `../python/pyalps/tests/fixtures/legacy_checkpoint.cpp`.
Frozen scheduler and solver results in `../python/pyalps/tests/fixtures/application_results/`
also exercise the public loaders. These prove historical ALPS compatibility;
they do not establish ALPSCore archive compatibility.

## Behaviors requiring a consolidation decision

The current suites characterize these historical behaviors:

- Forward-iterator list insertion into the fixed-capacity deque reverses the
  inserted list in the recorded scenario (`fixed_capacity_deque`, labeled
  `legacy_semantics`). This known defect must not constrain a replacement.
- Normalized autocorrelation is undefined (`NaN`) for constant samples. The old
  formatter substituted a display-only zero.
- `mcdata::mean()` includes the partial final bin; the free analysis function
  and `size()` operate on completed bins.

These are legacy characterization checks, not acceptance criteria for a new
statistics backend. Negation, nonnegative uncertainties and empty-histogram
construction instead have correctness assertions; known bugs are not contracts. Reconcile intended semantics explicitly, especially
correlated-sample error propagation. Add ALPSCore-produced archive fixtures and
parameter conversion checks, including integer widths, when evaluating ports.

The logarithmic `exp_number` tests express inequality as
`!(a == b)`: its missing inequality operator otherwise permits conversion of
distinct large values to indistinguishable infinities.

## Model, lattice, and Osiris semantic assertions

Basis tests enumerate ordered states, quantum numbers and fermion parity with
integer or half-integer references. Their fixture pairs remain distinct even
when they share an executable. Sign-problem cases preserve historical model
classifications; they do not prove arbitrary models sign-problem-free.

Lattice coloring checks the historical greedy assignment and verifies every
bond. The triangular lattice's four-color reference describes that ordering,
not an optimal coloring. Parity tests distinguish absent and empty
`BACKBONE_TYPES` and check site/bond metadata and periodic wrapping.

`osiris_dump_serialization` checks scalar widths, signs, large integers,
complex values and native/Boost adapter paths. Both historical readers consume
the immutable `xdrdump2.dump`; the writer must also reproduce its bytes,
independently of the current reader.

## Coverage boundaries

Timing programs and Alea `dumpbench.C` are benchmarks. The unregistered
`binned_data.C` refers to an absent header and provides no coverage. Remaining
manual MPI drivers and restart gaps are described in [the MPI guide](mpi.md).
Component guides record [HDF5 schema/type coverage](../src/alps/hdf5/tests/README.md),
[parameter contracts](../src/alps/params/tests/README.md) and
[graph coverage](../src/alps/graph/tests/README.md).
