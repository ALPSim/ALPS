# Retirement of the Boost.Python differential audit

The optional audit at revision `d5642c56fb95192510e03e228bc699f9286b87a0`,
under `script/pyalps_compatibility/`, compared 406 records, including API
inventories, against `f28d428017773f5794dc9896a544f95e8ef40443`. Its 19 allowed
differences were six kinds of intentional changes: full vector slices; integer
vector convergence arrays; signed-byte scalar/array storage; replacing a dataset
with an empty dictionary; extra max-num-binning result math methods; and removal
of Boost indexing implementation classes. It did not supply 406 independent
scientific correctness tests. Float normalization rounded to 11 significant
digits, inventories included binding implementation details, and probes recorded
exceptions as well as successful values. Its archived sources remain available
at that revision for a future targeted old/new investigation.

The audit is retired, rather than maintained as a second build of the obsolete
bindings. The following current tests cover its useful contracts (paths relative
to `python/pyalps/tests`):

| Audit area | Current evidence |
| --- | --- |
| Modules, public operations and mapping behavior | `test_binding_surface.py`, `test_mapping_lifetimes.py`; assert supported behavior, not every name from `dir()`. Boost indexing helper classes are intentionally absent. |
| Scalar/vector arithmetic, slices, copies and formatting | `mcdata_test.py`, `test_checkpoint_contracts.py`, `test_conversion_contracts.py`; historical values and correctness assertions replace rounded fingerprints. |
| Observable conversions and timeseries | `test_conversion_contracts.py`, `test_binding_surface.py`; integer convergence arrays, scalar/vector timeseries from archives, array layout acceptance and analytic timeseries statistics. |
| Signed-byte storage and rectangular/container conversions | `test_archive_dtypes.py`, `pyhdf5io_test.py`; dtype/shape assertions, exact large integers, empty-dictionary overwrite, readonly/Fortran/strided arrays. |
| Additional max-num-binning math | `accumulators_test.py::test_transcendental_functions_match_numpy` checks constant-bin analytic values for both mean and max-num-binning results. |
| Parameter values, ownership and persistence | `pyparams_test.py`, `test_binding_surface.py`, `test_checkpoint_contracts.py`; mutation, metadata, arrays and integer boundaries. |
| Historical checkpoint reads and RNG continuation | `test_checkpoint_contracts.py::test_historical_checkpoint_remains_readable` consumes immutable native-serializer bytes; `test_legacy_results.py` exercises historical application archives. |

This is deliberately not a claim of complete Boost.Python equivalence. The frozen
checkpoint does not exercise old Python object conversions, old bindings reading
new checkpoints, or the entire old seeded-RNG fingerprint corpus. Current RNG
round trips also cannot independently establish all historical seeded sequences.
Those are limits of the evidence, not promises silently supplied by the fixture.
If a specific compatibility requirement emerges, recover the relevant old probe
and baseline in matching interpreter/NumPy environments, then add a focused
regression. Do not restore inventories, environment-specific error strings, or
the complete legacy build just to increase test counts.
