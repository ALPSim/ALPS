# Accumulators

`ALPS::accumulators` is the accumulators library of the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore): `alps::accumulators`
with `<alps/accumulators.hpp>` and `<alps/accumulators/*.hpp>`. It is copied
from ALPSCore at commit `7146b9e1f017938a94e5dae35d88467cc5ba7969`
(2026-09-25).

It replaces ALPS's own NGS accumulators (`alps::accumulator`,
`<alps/ngs/accumulator/*.hpp>`) and the NGS observables over ALEA
(`alps::mcobservables`, `<alps/ngs/mcobservables.hpp>`). ALEA itself,
`<alps/alea.h>` in `alea/`, is unchanged.

## Building and linking

`ALPS::accumulators` is always built and is part of `ALPS::alps`. It needs
`ALPS::utilities`, `ALPS::hdf5` and the headers of `ALPS::numerics`, but not
BLAS or LAPACK:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS accumulators)
target_link_libraries(app PRIVATE ALPS::accumulators)
```

Accumulators and results are stored through ALPS's `alps::hdf5::archive`. In an
MPI build, `collective_merge` reduces through ALPSCore's
`<alps/utilities/mpi.hpp>` and `<alps/utilities/boost_mpi.hpp>` from
`ALPS::utilities`, and `ALPS::accumulators` links MPI.

The value types the accumulators are compiled for come from
`ALPS_ACCUMULATOR_VALUE_TYPES` in `<alps/config.hpp>`: `float`, `double`,
`long double` and `std::vector`s of them, as in ALPSCore.

## Provenance

The imported revision is under the MIT license in [LICENSE.TXT](LICENSE.TXT).
The imported copyright and license notices are retained in [COPYRIGHT.TXT](COPYRIGHT.TXT) and [LICENSE.TXT](LICENSE.TXT). The SDK installs both under `share/alps/licenses/accumulators/`. Publication pointers in imported comments have been removed; the source relationships below refer to the remaining code.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `accumulators/include/alps/accumulators.hpp` | `include/alps/accumulators.hpp` | Copied |
| `accumulators/include/alps/accumulators/**` | `include/alps/accumulators/**` | Copied, except that `convergence.hpp` and `feature/binning_analysis.hpp` gain includes (see below) |
| `accumulators/src/*.cpp`, `accumulators/src/feature/*.cpp` (11 files) | `src/`, `src/feature/` | Copied, except `mpi.cpp` (see below) |
| `accumulators/test/*.cpp` (31 files), `accumulators/test/accumulator_generator.hpp` | `tests/` | Copied, except two include lines (see below) |
| `utilities/test/vector_comparison_predicates.hpp` | `tests/vector_comparison_predicates.hpp` | Copied |
| `utilities/include/alps/testing/unique_file.hpp`, `utilities/src/unique_file.cpp` | `tests/support/` in the repository root | Copied; shared test support |
| `COPYRIGHT.TXT`, `LICENSE.TXT` | here | Copied |

Not copied:

- `accumulators/include/alps/accumulators_.hpp` and `accumulators/src/accumulators.cpp`,
  a "wrapped" interface ALPSCore installs but does not compile, and the test
  of it, `mean_wrapped.cpp`.
- The tests ALPSCore does not build: `divide_accumulators2.cpp`, and
  `acc_check_eigen.cpp`, which only checks that Eigen can be included.

The headers the accumulators need from ALPSCore's utilities module are in
`ALPS::numerics` and `ALPS::utilities`; see the
[numerics](../numerics/README.md) and [utilities](../utilities/README.md)
READMEs. The CMake files and this README are ALPS's own.

### Local changes

Two imported headers relied on their includers for some declarations. Each now
includes what it needs, so it compiles on its own, and a comment marks the
addition:

- `alps/accumulators/convergence.hpp` includes `<boost/throw_exception.hpp>`
  and `<stdexcept>`.
- `alps/accumulators/feature/binning_analysis.hpp` includes
  `feature/error.hpp`, whose `error_type` it uses.

One source change, for MPI builds:

- `src/mpi.cpp` dispatches its reduction helpers on
  `std::integral_constant<bool, hdf5::is_content_continuous<T>::value>`, and
  likewise for `is_continuous`, instead of on the traits' own `::type`. ALPS's
  `alps::hdf5` traits derive from `boost::true_type` and `boost::false_type`,
  ALPSCore's from the `std::` types the helpers take.

One change to `FullBinningAccumulator` results, marked `ALPS:` in
`feature/max_num_binning.hpp` and `src/feature/max_num_binning.cpp`:

- Adding, subtracting, multiplying or dividing by a constant, or negating,
  shifts or scales the result's mean and error. `-r`, `r + 1` and `r / 2` have
  the error of `r`, and `2 * r` twice it. ALPSCore replaced the
  binning-analysis error with a jackknife over the bins, which for correlated
  data differs from it: by 7.5% for an AR(1) process with coefficient 0.8.
  Operations on two results and nonlinear functions still use the jackknife.
  `tests/linear_transform_error.cpp`, ALPS's own, tests this.

One test change:

- `tests/save_load.cpp` and `tests/mean_err_count.cpp` include
  `vector_comparison_predicates.hpp` from their own directory. ALPSCore names it
  by its relative path in ALPSCore's source tree,
  `../../utilities/test/vector_comparison_predicates.hpp`.

The test registration also differs:

- ALPS runs each GoogleTest case as its own CTest test, so the `save_load2`
  cases, which write the same files, share a CTest resource lock.
- The MPI tests (`mpi_merge`, `mpi_merge_uneven`, `repeated_merge` and
  `zero_vector_mpi`) keep their own `main()` and run on 2 and 3 ranks. ALPSCore
  ran them on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- The tests' `TYPED_TEST_CASE` deprecation warnings are silenced.

`tests/ar1_calibration.cpp` is ALPS's own. It checks the errors and
autocorrelation times against the exact values of an AR(1) process.

## Known limitations

These are properties of the imported ALPSCore code; ALEA's observables share
the first.

- The variance is accumulated as the sum of squares minus the squared sum.
  It loses precision when the mean is large compared with the spread: for
  unit-variance noise on an offset of 10⁶, the log-binning error is about 25%
  too small, and from about 10⁷ it is NaN. Constant data often gives a NaN
  error for the same reason.
- `LogBinningAccumulator` reports an infinite error below 256 samples.
  `FullBinningAccumulator` reports the unbinned error there.
- Arithmetic on `NoBinningAccumulator` and `LogBinningAccumulator` results
  propagates errors as if the operands were uncorrelated. For correlated
  observables, such as ⟨m⁴⟩/⟨m²⟩², use `FullBinningAccumulator`, whose results
  use the jackknife.
- `float` accumulators lose the mean and error after about 10⁶ samples.
- `FullBinningAccumulator`s cannot be merged.
