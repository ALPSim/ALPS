# Reproducing the Python migration evidence

These checks use native Python environments and the existing ALPS SDK. They
do not need containers, a virtual machine, or a new CI build matrix. The normal
wheel smoke jobs run the same binding suite as before and upload its evidence.
Legacy comparisons and the six solver workflows are opt-in local checks.

## Installed-wheel validation

Install a wheel in an environment containing NumPy, SciPy, and pytest. From
the repository root run:

```sh
python script/validate_pyalps.py --output _build/validation --wheelhouse wheelhouse
```

The output contains a JUnit report, test log, source revision and dirty-file
hashes, installed package/dependency versions, extension hashes, distribution
hashes, executed commands, timings, and pass/fail status. A failing check exits
nonzero and retains its evidence. Each subprocess has a five-minute timeout.
Outputs belong in an ignored build directory or outside the checkout.

Use `--packaging` to include release/packaging tests (requires Python 3.11 or
newer, packaging, and scikit-build-core). Use `--downstream` to enable the two compiled consumers;
this requires the matching SDK installed at `_build/wheel-deps/install`, a
C++ compiler, CMake, and nanobind. Check the reported skips: a standard wheel
smoke run does not exercise these opt-in consumers or MPI without mpi4py.

Use `--applications` for input generation, dispatch, and result loading through
the installed `spinmc`, `loop`, `dirloop_sse`, `sparsediag`, `fulldiag`, and
`dmrg` programs. The exact-diagonalization cases assert the four-site periodic
Heisenberg ground-state energy, -2, within 1e-10. The Monte Carlo and short
DMRG cases assert finite loaded results; they do not establish convergence.
Each run writes its own inputs and outputs under the evidence directory.

## Old/new behavioral comparison

The reference is master at `f28d428017773f5794dc9896a544f95e8ef40443`.
Build its Boost.Python modules in a separate checkout with the **same Python
and NumPy versions** as the new-wheel environment. NumPy 1.26 and Python 3.10
are a useful common baseline. Keep old and new modules in separate processes.

`legacy-build.patch` records the build-only adaptations used in the original
audit: expose the experimental accumulator module without selecting the
unrelated global new-ALEA mode, and link Python modules through `Python::Module`
without standalone Python's redundant extra libraries. It does not patch the
compared binding algorithms. Apply it only in the legacy checkout.

For example, with `LEGACY_SOURCE`, `LEGACY_BUILD`, `LEGACY_PYTHON`, and
`PR_SOURCE` set to your paths:

```sh
git -C "$LEGACY_SOURCE" checkout f28d428017773f5794dc9896a544f95e8ef40443
git -C "$LEGACY_SOURCE" apply "$PR_SOURCE/script/pyalps_compatibility/legacy-build.patch"
cmake -S "$LEGACY_SOURCE" -B "$LEGACY_BUILD" \
  -DALPS_BUILD_PYTHON=ON -DALPS_BUILD_APPLICATIONS=OFF \
  -DALPS_BUILD_TESTS=OFF -DALPS_BUILD_EXAMPLES=OFF \
  -DALPS_BUILD_DEVELOPER_TOOLS=OFF -DALPS_BUILD_FORTRAN=OFF \
  -DALPS_ENABLE_MPI=OFF -DALPS_ENABLE_OPENMP=OFF \
  -DALPS_NGS_USE_NEW_ALEA=OFF -DPython_EXECUTABLE="$LEGACY_PYTHON"
cmake --build "$LEGACY_BUILD" --parallel 2
```

Use your normal HDF5 and Boost configuration. The original audit used Boost
1.87, Apple Clang, and a package-specific SZIP include directory to avoid
mixing Homebrew Boost headers with the selected Boost sources. This baseline
build is a one-time local prerequisite, not an additional CI job.

Then run from the PR checkout, using the new-wheel environment's Python:

```sh
python script/validate_pyalps.py --output _build/compatibility \
  --legacy-python "$LEGACY_PYTHON" \
  --legacy-modules "$LEGACY_BUILD/lib/pyalps"
```

`probe.py` records public names, class members, scalar/vector arithmetic,
observables, timeseries, RNG sequences, HDF5 types/shapes/layouts, and parameter
values. It includes the three rectangular-table cases found in the adversarial
review. The original corpus had 403 records; this corpus has 406. Records
include inventories and probes, not 406 independent numerical tests.

`compare.py` rejects mismatched Python/NumPy versions and unexpected changes.
`expected_differences.json` contains the exact old and new values and reasons
for the 19 previously examined differences. It is not a name-only allowlist:
a different value under an allowed name still fails. A disappearing expected
difference also fails and must be examined before updating the file.

The numeric fingerprints round floating values to 11 significant digits, as
in the original corpus, to suppress insignificant floating-point noise. This
is not a universal scientific tolerance or a substitute for solver-specific
reference tests. Arrays retain dtype and shape; selected parameter probes also
record their Python type. The ordinary tests assert arithmetic and ownership
contracts that a normalized fingerprint alone cannot establish.

`test/pyalps/test_archive_dtypes.py` also checks mixed numeric rows at the
float64 precision boundary and the int64/uint64 limits, in both row orders and
with nested containers. These tests compare integers as Python ints: ordinary
NumPy int/float equality can itself round the expected integer and hide data
loss. Lossy mixtures must retain separately typed rows; exactly representable
mixtures must still load as arrays. These fast cases run in the existing wheel
smoke suite and need no additional CI job.

The runner also writes checkpoints using each implementation and reads them
with the other. These cover ordinary parameters, RNG continuation, scalar and
vector observations, and persisted results. Legacy `result.load()` is broken;
the legacy side checks the stored result datasets, while the new side performs
actual result restoration. This does not claim compatibility for arbitrary
metadata or all native parameter dtypes.

## Local sdist rebuild and additional checks

Build the matching native SDK using the documented `wheel-deps` preset. Set
`ALPS_DIR` to its installed `share/alps` directory. Then, outside the checkout:

```sh
python -m build --sdist --outdir dist "$PR_SOURCE/bindings/python/pyalps"
mkdir extracted
tar -xzf dist/pyalps-*.tar.gz -C extracted
python -m build --wheel --outdir wheelhouse extracted/pyalps-*
```

Repair the wheel with the normal platform tool (delocate on macOS, auditwheel
on Linux), install it in a fresh environment, and run `validate_pyalps.py`
against that environment with `--wheelhouse` pointing at the repaired wheel.
This is deliberately a local release check: the normal CI sdist job checks
metadata and payload, while the existing wheel jobs already perform full builds.

MPI remains covered by the existing CI job. A local two-rank run is:

```sh
mpiexec -n 2 python -m pytest -q \
  test/pyalps/test_binding_surface.py::test_mpi4py_compatibility_surface \
  test/pyalps/test_mpi_requests.py
```

The native-only parameter loader still has an inherited limitation: unsupported
numeric storage types such as a standalone int64 or float32 dataset can leave
the default value at zero. Python-owned parameter loading uses a different
decoder. The migration evidence must not be read as a blanket claim that every
native checkpoint dtype is handled. This issue is separate from the repaired
rectangular-table regression.
