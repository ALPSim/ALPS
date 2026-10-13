# HDF5 test contracts

All registered runtime tests use GoogleTest and are discovered by CTest under
`hdf5.<target>.*`. Each case owns a temporary directory, so parallel runs do not
share archives. Serialization comparisons are exact: they check lossless storage
of a value already represented in memory, not a numerical approximation.

## Migration inventory

| Legacy executable/scenario | GoogleTest contract |
| --- | --- |
| `hdf5_complex`, `hdf5_bool` | Scalar/vector complex values and nested bool, enum, pair, unsigned values survive storage. Previously these largely printed values. |
| `hdf5_copy` | Recursive copy preserves ragged vectors, complex data, attributes, and caller contexts; unrelated source datasets are not copied. |
| `hdf5_real_complex_vec`, `hdf5_real_complex_matrix` | Reading real storage into complex containers rejects missing complex metadata. This is the current explicit contract in `vector.hpp` and `matrix.hpp`; the former tests incorrectly accepted both conversion and rejection. |
| unregistered `hdf5_real_complex.cpp` | Removed redundant combined source; both scenarios are registered independently in the preceding two tests. |
| `hdf5_replace`, `hdf5_pair` | Writable archive mode retains existing datasets; real pointer storage can be overwritten by complex pointer storage and read back. |
| `hdf5_family`, `hdf5_multiarchive` | Multiple simultaneously open handles and subsequently reopened handles read the same stored value. |
| `hdf5_valgrind` | Repeated archive open/write/read/close exercises lifetime management and checks each value; sanitizer jobs can detect memory defects. |
| `hdf5_multi_array`, `hdf5_vecvecdbl`, `hdf5_vecveccplx` | Uniform/ragged shapes and all stored elements are preserved. |
| `hdf5_memory` | Complex-to-real read is rejected; the same archive remains usable for a subsequent real read. |
| `hdf5_misc`, `ngs_hdf5` | Scalar conversions, vectors, initialized raw arrays, user-defined serialization, strings, attributes, metadata and overwrites retain their values. |
| `hdf5_exceptions` | Missing dataset throws and its first diagnostic line names the path; the old `.output` expectation is now an assertion. |
| `hdf5_fortran_string` | Historical binary fixture remains unchanged; both its string datatype and the independently inspected value `N_total` are checked, including the legacy reader's trailing NUL in its returned `std::string`. |
| `hdf5_omp` | Each OpenMP worker writes and reads its own archive. Exceptions are captured within workers and reported in the main test thread. |
| `type_check.cpp.in` | The existing type list and its dataset/attribute/SZIP modes remain. Round trips and overwrite scenarios are separate discovered cases, each retaining 32 iterations. |
| formerly unregistered `hdf5_large` | An explicit `ALPS_BUILD_LARGE_MEMORY_TESTS=ON` stress test preserves the large family-file writes and checks dimensions/endpoints. It requires a 2 GiB allocation and approximately 4 GiB of disk and is off by default, independently of extensive tests. |

The generated type matrix reports unavailable SZIP encoding and unsupported
attribute combinations as named GoogleTest skips. Shared-array overwrite had no
legacy implementation and is not registered; shared-array round trips
still cover one-, two-, and three-dimensional storage. Scalar overwrite retains
all 32 iterations and concurrent handles from the actual legacy execution. Nine
unreachable cross-type transitions following its first `return` were removed.
They cannot be enabled indiscriminately: attributes do not resize their storage,
and group serializers do not replace an existing scalar dataset. Supported
real/complex replacement remains covered by `hdf5_pair` and `ngs_hdf5`.
Random/empty/special values and pointer
shape comparisons report their own failures instead of collapsing to an exit
status. The deterministic generator retains seed 42.

Cross-component tests in `tests/integration/hdf5` check legacy parameter type
changes, observable-set count/mean restoration, and all eleven Ising inverse
temperatures from 0 to 1. Ising checks compare archived values and parameters to
the simulation, verify finite uncertainties and exact physical bounds, and do
not introduce probabilistic pass/fail thresholds.

## ALPSCore conformance

Tests ported from ALPSCore 7146b9e1 keep ALPSCore's copyright header and name
their source and every change in a comment at the top. ALPSCore is MIT
licensed; its copyright lines are
"Copyright (C) 2012-18 University of Michigan" and "Copyright (C) 1998-14 ETH
Zurich", and it asks publications to cite Wallerberger et al.,
arXiv:1811.08331.

| Test | Source | ALPS changes |
| --- | --- | --- |
| `hdf5_attributes` | `hdf5/test/hdf5_attributes.cpp` (with fix c7c64ca0) | ALPS temporary directory; cases for ALPS's `path/@name` spelling and a root-context attribute |
| `hdf5_open_strmode` | `hdf5/test/hdf5_open_strmode.cpp` (with fix 33e26902) | ALPS temporary directory; `"wl"` is valid in ALPS, so `"wx"` and `"R"` are the invalid modes; a case for combined letters |
| `hdf5_native_type` | `hdf5/test/hdf5_native_type.cpp` (with 4daa2da2) | none |
| `hdf5_boost_optional` | `hdf5/test/hdf5_boost_optional.cpp` (with 7b5b7be5) | ALPS temporary directory; fixed values instead of ALPSCore's `testing::datapoint`; `TYPED_TEST_SUITE`; a case reading into an empty optional |
