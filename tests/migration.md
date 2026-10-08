# C++ test migration inventory

## Numerics, containers, and alea

This inventory maps the registered tests on the PR #163 base (`fee875ce`) to
GoogleTest. The work is carried forward onto the updated PR head `25702d9`.
CTest discovers each runtime case. The three components have 689
cases: 647 numeric, 17 fixed-capacity, and 25 alea. Counts describe execution,
not the amount of behavioral coverage: some cases contain a long sequence of
related operations.

## Preserved scenarios

| Previous executable(s) | GoogleTest coverage |
| --- | --- |
| `matrix_unit_tests` | All 30 scalar scenarios across float, double, int, unsigned int, unsigned long, complex float, and complex double; mixed matrix/vector multiplication across all four original type pairs. Construction, iterators, resizing, allocation-failure guarantees, insertion/removal, arithmetic, multiplication, and conjugation/transposition remain covered. |
| `matrix_column_view`, `matrix_transpose_view`, `matrix_ublas_sparse_functions` | All 11 column-view, nine transpose-view, and one sparse-multiplication scenarios across the same seven types. |
| `matrix_hdf5` | Scalar and nested-matrix serialization across seven types, reserved-capacity behavior, column-major disk axes, complex vector layout, and empty-matrix loading. Each case owns a temporary directory. |
| `matrix_deprecated_hdf5_format_test` | The original immutable HDF5 input, all 60 matrix entries, dimensions, and capacity. Expected decimals come from the former output fixture; absolute bounds match its six-significant-digit precision. |
| `matrix_algorithms` | All 13 algorithms/invariants across all 18 original shape/type combinations: trace, transpose, adjoint, SVD, LQ, QR, and inverse. The original Boost percentage tolerance of `1e-6` becomes a relative bound of `1e-8`; identity residuals retain the `1e-6` absolute bound. RNG seed 3 is reset for each case so selection and execution order cannot change inputs. |
| `matrix_kron` | The four original equal/unequal rectangular operand shapes, individually parameterized, checking every output entry and both dimensions. |
| `real_tests` | Qualified and unqualified lookup for all 12 scalar, complex, vector, and nested-vector types; recursive independent real-part references replace printed values. |
| `vector_functions` | Every original arithmetic and transcendental expression, compared against scalar standard-library references for all ten elements. |
| `vector_valarray_conversion`, `accumulate_if` | Both conversion directions, the element-type conversion overload, and both accumulation overloads. Empty/rejected ranges are additional boundary checks; arithmetic-series identities supply accumulation references. |
| `fixed_capacity_vector`, `fixed_capacity_deque` | The original modifier, range-insertion, iterator, copy, assignment, and different-size swap scenarios, now asserting container contents directly. |
| `fixed_capacity_traits` | All six dynamic-capacity and five fixed-capacity types, including container adaptors. |
| `test_vector`, `test_deque` | The original 10,000-operation traces against standard containers, with explicit historical seed 331. Both plain and non-POD values are checked after each operation. External live-address tracking replaces undefined reads of uninitialized fields in the old lifetime test double. |
| `detailedbinning`, `simpleobseval`, `vectorobseval`, `mcdata` | Original measurement order, sample counts, seed 1, thermalization, construction/assignment paths and evaluator names, merges, slicing, and transformed expressions. Every numerical estimate, printed bin-error level/count, and defined autocorrelation value has a direct assertion. |
| `mcdata2` | Every active scalar and vector arithmetic/assignment/transcendental expression, in separate cases. Means and errors retain the original ten-significant-digit reference precision, including every element previously abbreviated in vector output. |
| `mergeobs_*` | All six merge/transform scenarios, each with two sets of 128 samples and seed 1. The old `.op` reference values now assert means, errors, counts, convergence states, and available autocorrelation values. These files previously were not selected by the `.output` runner. |
| `complexobservable` | The original 1,000-sample complex-observable/product scenario, now checking independent sums and finite product statistics. An analytic constant-product/zero-error case is additional coverage. |
| `mcdata_transform_variance` | All four existing regression scenarios, independently discoverable: supplied variance for linear, unary, and binary transforms, and clearing linear-transform variance. |
| `observableset_hdf5` | Saving/loading empty sets, populating the loaded set, ratio construction, and loading into both empty and prepopulated sets. All prior estimates, histogram counts, detailed error levels, and defined tau values are asserted; unexpected exceptions now fail. |
| `mcanalyze` | The original archive/analysis workflow, with independent all-sample and completed-bin means. `mcdata::mean()` includes the partial final bin; the free analysis function and `size()` operate on 78 completed bins of 128 samples. |
| `histogram2` | Merging the two empty histograms, followed by direct checks of the merged range, zero count, and zero-filled bins. |
| `histogram`, `signed`, `testobservableset`, `observableset_xml` | Real GoogleTest cases retain the exact XML fixtures and original inputs through `StreamFixture`. Formatting/schema compatibility remains intentional here. |

The alea numerical reference tables are transcribed from the base commit's
fixtures, not regenerated from the implementation under test. Their tolerances
express historical printing precision; they are deterministic regression checks,
not statistical acceptance intervals. Obsolete numeric stdout fixtures are
removed. XML fixtures and the historical matrix HDF5 input remain versioned.

Sanitizer validation exposed two pre-existing production defects in these
components. Empty vector/valarray conversions accessed element zero before
performing an empty copy; all three container overloads now avoid that access.
`ObservableSet::clear()` erased owning pointers without deleting their objects;
it now applies the same deletion operation as the destructor before clearing the
map and sign metadata. The existing signed/unsigned XDR reload scenarios exercise
that lifetime path with leak detection enabled.

The same run found unreclaimed scalar HDF5 variable-length strings in the archive
reader. Dataset and attribute reads now own that HDF5 allocation through
conversion, including exceptional conversion paths. These ownership and boundary
fixes are separate from the numerical compatibility decisions below.

## Behaviors requiring a consolidation decision

The migration preserves these observed historical behaviors without changing
production code:

- Forward-iterator list insertion into the fixed-capacity deque reverses the
  inserted list in the recorded scenario.
- Unary negation of `mcdata` retains the positive mean in the original scalar
  and vector scenarios.
- Multiplying/dividing vector evaluators by negative constants can produce
  negative reported errors.
- Normalized autocorrelation is undefined (`NaN`) for constant samples. The old
  text formatter substituted a display-only zero; direct tests distinguish that
  presentation behavior from the numerical API.

An additional empty-histogram evaluator construction exposed a crash while
expanding the empty-merge assertions. That conversion was not exercised by the
original test. The migrated case inspects the already-merged histogram directly;
the additional constructor path needs a separate production regression/fix.
These are review topics for the compatibility contract, not endorsements of the
current numerical semantics.

## Unregistered sources

- `containers/tests/timing_{queue,stack,vector}.C` and alea `dumpbench.C` are
  benchmarks, not correctness tests; this migration does not claim their coverage.
- Alea `binned_data.C` was unregistered and references the absent
  `alps/alea/binned_data.h`. It remains an obsolete source, not a skipped passing
  test.
- Alea `observableset_mpi.C` is now registered at two and three ranks; its
  executed transport contracts are described in [native MPI tests](mpi.md).

## Validation

Build the component test targets, then run:

```sh
ctest --test-dir <build> -L '^(numeric|fixed_capacity|alea)$' \
  --parallel 4 --output-on-failure --output-junit component-tests.xml
```

All seven numeric scalar types, all four mixed-type pairs, and all 18 algorithm
shape/type combinations must appear in discovery. A successful executable build
alone does not establish that these cases ran.

## Legacy parameters and parapack

The registered serial tests in these components now use GoogleTest: 45 cases in
nine legacy-parameter executables and 35 cases in 16 parapack executables. Focused
runs passed all 80 cases. The formerly unregistered `parameters_hdf5` and
`id2string` programs contribute one and nine of those cases, respectively.

| Previous executable(s) | Preserved or activated coverage |
| --- | --- |
| `expression`, `expression2`, `flatten` | Historical parameterized evaluations, unresolved references, seeded random-function results, canonical simplification/flattening strings, and numerical equivalence. Known analytic values supplement the old decimal references. |
| `parameter`, `parameters`, `parameterlist` | Quoting, whitespace, environment substitution, key/value ordering, replacement/deletion/copying, inherited parameter blocks, `#clear`/`#stop`, and an isolated XDR round trip. The previous inputs are embedded in the semantic cases. |
| `parameters_xml`, `parameterlist_xml` | The original input/output fixtures remain exact XML spelling and ordering contracts through `StreamFixture`. |
| `parameters_hdf5` | Newly registered HDF5 round trip of the historical numbers, paths, trailing whitespace, and unevaluable strings, using a unique temporary directory. |
| `clone_info`, `clone_phase` | Existing exact serialized XML transcripts and XDR/HDF5 round trips, with isolated files. Their fixtures remain versioned. |
| `clone_timer`, `info_test`, `time` | Explicit clock values for timer/priority behavior; clone identity, seeds, phases, progress, checkpoint identity, and bounded timestamp checks replace busy loops and sleeps. |
| `exmc_optimize`, `exp_number`, `linear_regression`, `merge`, `wl_weight` | Historical ladders/interpolation, signed logarithmic arithmetic, independently centered regression references, merged means/errors, and seeded Wang–Landau visit/weight results. Numeric expectations are direct assertions rather than formatted stdout. |
| `footprint`, `integer_range`, `percentage`, `temperature_scan`, `version` | Object/container accounting, all four original integer domains and parsing/scaling/range operations, percentage spellings, the original 47-step temperature/progress schedule, and release-version syntax. |
| `id2string` | All nine historical process identifier spellings are individually discovered string assertions. |

Removed numeric/text fixtures are replaced by explicit values or independent
calculations in these cases. The exact XML input/output fixtures remain;
obsolete empty MPI inputs and numeric MPI transcripts are removed after their
assertions are activated. The unused `expression.input2` exploratory complex
input remains outside registered coverage.

Numeric tolerances preserve the original printed precision where that is the
available reference, rather than treating rounded decimals as exact values.
The merged mean's `1e-8` absolute bound covers accumulation roundoff and is
stricter than the former three-significant-digit output. These are compatibility
checks, not independent validation of statistical uncertainty.

The historical logarithmic `exp_number` type has an equality operator but no
inequality operator. The migrated tests express inequality as `!(a == b)`;
otherwise implicit conversion to `double` can turn two distinct large values
into infinities before comparison. The migration preserves the logarithmic
comparison contract without changing production operators.

See [native MPI tests](mpi.md) for newly activated parameter/observable transport,
parapack process/locking/metadata tests, scheduler tests, and the explicitly
unregistered manual Ising checkpoint/restart drivers.

## Model, lattice, and Osiris semantic assertions

The follow-up to `a357f73ad` replaces 13 transcript cases with 38 semantic cases
in the same executables. Their expected values come from the pre-existing
fixtures at that commit; none were regenerated from current test output.
The removed input files are represented explicitly by named cases and parameters.

| Executable(s) | Preserved contracts and additional checks |
| --- | --- |
| `model_example2` | Six cases cover fermion occupations, hardcore bosons, and spin 1/2, 1, 3/2 and 2. Quantum-number names, ordered states, and exact half-integer values replace printed basis listings. Index lookup and sortedness are also checked. |
| `model_example6` | Four cases retain default and spin-2 spinful bosons (occupation cutoff 2), t-J and alternative t-J. The spinful reference enumerates integer `(N, J, Jz)` tuples with `0 <= J <= N*spin` and `-J <= Jz <= J`; the explicit 5- and 35-state contracts are unchanged. |
| `model_example7` | Separate 14-state and 55-state cases preserve enumeration before and after quantum-number descriptor addition. The combined range reaches occupation 4; this is not a tensor-product state-count claim. |
| `model_example8` | All eight historical sign-problem classifications are separately discoverable with named lattice, exchange and field parameters. Unspecified diagonal coupling and field still use model defaults. These are historical classifications, not a general proof that arbitrary models have no sign problem. |
| `model_example9` | Hardcore boson and spinless fermion single-quantum-number representations retain their ordered empty/occupied states. Fermion parity is checked explicitly too. |
| `coloring` | All four lattices retain the exact historical greedy color assignment and count. Every bond must also connect different colors, and all color indices must be valid. The four-color triangular result describes this vertex ordering, not an optimal coloring. |
| `parity` | All six backbone scenarios retain site parity, including the difference between absent and empty `BACKBONE_TYPES`. Direct assertions preserve graph dimensions, site types/coordinates, bond indices/endpoints/types/displacements and periodic wrapping previously present in the XML. Bipartiteness and the public numeric parity accessor are checked too. |
| `xdrdump`, `xdrdump2`, `boostdump`, `boostdump2`, `boostdump3`, `boostdump4` | Six cases retain scalar width/sign, both large signed-integer scenarios, strings, complex values, and native/Boost adapter paths. Double values now compare at full precision instead of six printed significant digits. Both historical readers still consume the immutable `xdrdump2.dump`; adapter output must also match its exact bytes, independently of the current reader. |

Private helpers in `model/tests/site_basis.cpp` report failures by state and
quantum-number name; they do not extend the shared testing framework. Numerical
enumeration references use integer/half-integer arithmetic and need no tolerance.
Exact XML serialization and remaining symbolic-expression fixtures are retained
where their representation is the compatibility contract. Those cases still need
additional semantic coverage when their underlying behavior is changed.

## Test suite consolidation

Related cases now share a source file and executable where they have the same
dependencies and execution requirements. Every GoogleTest suite/case name,
assertion, numerical tolerance, and reference fixture is retained. CTest still
discovers and executes cases separately; only the target portion of the CTest
name changes for the following groups:

| Previous targets | Current target | Cases |
| --- | --- | --- |
| `model_example2`, `model_example6`, `model_example7`, `model_example9` | `model_site_basis` | 14 |
| `model_example10`, `model_example12` | `model_basis_states` | 2 |
| `model_example15`, `model_example16` | `model_bloch_basis_states` | 2 |
| `xdrdump`, `xdrdump2`, `boostdump`, `boostdump2`, `boostdump3`, `boostdump4` | `osiris_dump_serialization` | 6 |

The basis-state and Bloch-state pairs shared identical code but different
inputs and expected outputs. Both fixture pairs remain unchanged and each is
checked by its own case. The dump suite shares record writing and value checks
while preserving native conversion/getter APIs, Boost adapter APIs, extended
integer/complex values, historical-file readers and the byte-for-byte writer check.

Repeated catch-and-fail wrappers are removed from 28 model, lattice and XML
serialization cases. GoogleTest reports uncaught exceptions as failures;
`StreamFixture` still restores redirected streams during unwinding. Inner scopes
still destroy serialization objects before comparing output. This cleanup changes
no fixture contents and adds no shared framework machinery.

Validation retained all 918 native CTest cases after applying the target mapping.
All 24 consolidated cases also passed two shuffled in-process iterations to check
for shared-state and order-dependent failures; recorded fixture hashes were unchanged.
