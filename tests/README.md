# Testing ALPS

C++ tests use GoogleTest and run through CTest; Python tests use pytest. Component
tests live beside their code in `src/alps/<component>/tests`, Python tests in
`python/pyalps/tests`, and repository-wide integration, installed-library,
packaging and CI checks here. The two suites are independent: passing CTest says
nothing about the Python package.

## When to add a test

A test earns its place when a real bug would make it fail. A bug fix adds the
one test that fails without the fix; a new feature checks its result against an
exact or independent value; a refactor adds none. A test that would only fail
because the code changed, not because a result is wrong (a default value read
back, one case per overload), is removed in review.

## Local workflows

With the prerequisites from [CONTRIBUTING.md](../CONTRIBUTING.md#prerequisites):

```sh
cmake --preset dev
cmake --build --preset dev --parallel 2
ctest --preset dev                 # whole native suite
ctest --preset dev -L hdf5         # by label
ctest --preset dev -R Matrix       # by name
```

Presets: `dev` is the library only; `default` adds the applications; `mpi`
enables MPI (`ctest --preset mpi -L '^mpi$'` runs the rank topologies);
`debug` builds applications and tests with assertions on; `sanitizers` adds
ASan and UBSan; `extensive` enables the graph and HDF5 type-matrix tests. Build
a preset before testing it. GoogleTest is a test-only dependency, found if
installed (>= 1.14) or fetched (pinned 1.18.0; offline:
`-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/path`).

Native MaxEnt reference tests need NumPy and h5py in the CMake-selected
Python. Python suites, after the editable install in CONTRIBUTING.md:

```sh
python -m pytest python/pyalps/tests -q -rs
python -m pytest tests/cmake -q -rs                            # installed-library consumers; needs ALPS_DIR
ALPS_XML_BUILD="$PWD/_build/default" python -m pytest tests/cli -q   # XML CLI; needs xsltproc
python -m pytest tests/ci tests/packaging -q                   # workflow and release helpers; Python >= 3.11
python -m pytest tests/tutorials -q                            # standalone numerical tutorials; numpy only
```

`PYALPS_TEST_DOWNSTREAM_EXPORT=1` adds the slower tests that compile Python
extensions against the installed library; `tests/cmake` expects the MPI-disabled
LP64 build of the `default` preset, and `ALPS_TEST_CMAKE_ARGS` passes a JSON
array of toolchain arguments to its temporary projects. Read the skip reasons
(`-rs`): a skipped case is not coverage. `python
.github/scripts/validate_pyalps.py --output _build/validation --wheelhouse
wheelhouse` records installed-wheel validation evidence (`--packaging`,
`--downstream`, `--applications` add checks). Keep
`python/pyalps/tests/fixtures/legacy_checkpoint.h5` frozen: regenerating it
with the current library would lose the backward-compatibility check.

## Writing a C++ test

Make an ordinary executable that links only the component under test and
register it with `alps_add_gtest(target COMPONENT hdf5 LABELS compatibility)`
in the test `CMakeLists.txt`; every case becomes a CTest entry. Group cases in
one executable when they share dependencies and labels; keep separate targets
for MPI, stress or compiler-specific tests. Labels: `compatibility`,
`integration`, `scientific`, `extensive`, `mpi`, and `legacy_semantics` for
characterisation of old behaviour that a replacement need not keep. `TIMEOUT`
(default 60 s, MPI 120 s) and `THREADS` are set per target when justified.

Cases in one executable may run concurrently, so use
`alps::testing::TemporaryDirectory` for files and never depend on ordering or
another case's output. Use GoogleTest assertions, not `assert()`, which Release
builds compile out. Textual contracts use `StreamFixture`
(`tests/support`) and compare inside the test, not by diffing printed output.

MPI tests use `alps_add_mpi_gtest` and run as complete rank topologies; every
rank must reach the same collectives even after a non-fatal assertion failure.
See [mpi.md](mpi.md).

## Choosing assertions and reference data

- Assert values directly. For floating point, state the tolerance and where it
  comes from; do not widen a tolerance to pass a new platform. Tolerance
  assertions reject NaN, so test intended non-finite behaviour explicitly.
- Prefer analytic results, independent calculations and invariants
  (normalisation, conservation) over values the code itself produced. A
  writer/reader round trip can hide a shared bug; keep historical files too.
- Seed stochastic tests and print the seed. Such a test reproduces one
  trajectory; it does not validate error bars. New statistical validation
  states its ensemble, autocorrelation treatment and false-failure probability.
  Do not re-run a failed numerical test until it passes.
- Fixture changes carry provenance and a behavioural explanation in the PR.
- Skip with a named reason when an optional capability is missing. Do not add
  permanently disabled tests; a quarantine names the reason, issue and owner.

Historical references migrated from the Boost.Test suites, the quirks kept
pending a consolidation decision, and the per-component inventories are in
[migration.md](migration.md), [mpi.md](mpi.md),
[hdf5](../src/alps/hdf5/tests/README.md) and
[graph](../src/alps/graph/tests/README.md).

## CI

Pull requests report one aggregate `CI` check that selects the affected areas
(see [CI coverage](../CONTRIBUTING.md#ci-coverage)); weekly and release runs
cover the compatibility matrix. Timing programs are not correctness gates and
shared-runner wall times are not regression thresholds.
