# ALEA v2

`ALPS::aleav2`, ALEA v2, is the ALEA library of the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore): accumulators for
mean, variance, covariance, autocorrelation and batch statistics, error
propagation, and result serialization. It is copied from ALPSCore at commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25). Its include spellings
are ALPSCore's: `<alps/alea.hpp>` and `<alps/alea/*.hpp>`, in namespace
`alps::alea`.

ALPS's own ALEA, in [alea/](../alea/) with `<alps/alea.h>` and
`<alps/alea/*.h>`, remains the default. It is unchanged and not deprecated.
Nothing in ALPS uses ALEA v2 yet, and the two will be swapped only after a
conformance suite passes on both. Only the names ALPS owns say v2: this
directory, the target `ALPS::aleav2`, the SDK component `aleav2` and the test
names. The two libraries share the include prefix and namespace without
sharing a file name, and one program can use both; `aleav2_coexistence`
checks this, including both kinds of result in one archive.

## Building and linking

ALEA v2 needs [Eigen](https://eigen.tuxfamily.org) 3.3 or newer; it was tested with
Eigen 5.0.1. ALPS builds `ALPS::aleav2` only when CMake finds Eigen, and reports
otherwise during configuration. The applications need Eigen for MaxEnt, so only
a library-only build (`ALPS_BUILD_APPLICATIONS=OFF`) can leave ALEA v2 out, with
`CMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON`. An installed SDK provides it as a
separate component, so other consumers do not need Eigen:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS aleav2)
target_link_libraries(app PRIVATE ALPS::aleav2)
```

ALEA v2 stores results through ALPS's `alps::hdf5::archive`
([`hdf5_serializer`](include/alps/alea/hdf5.hpp)).

In an MPI build, `<alps/alea.hpp>` also includes ALPSCore's
[`mpi_reducer`](include/alps/alea/mpi.hpp), which combines accumulators across
ranks through `<alps/alea/mpi.hpp>` and ALPSCore's `<alps/utilities/mpi.hpp>`
from `ALPS::utilities`; `ALPS::aleav2` then links MPI. Serial builds install
neither header.

## Provenance

The imported revision is under the MIT license in [LICENSE.TXT](LICENSE.TXT).
The imported copyright and license notices are retained in [COPYRIGHT.TXT](COPYRIGHT.TXT) and [LICENSE.TXT](LICENSE.TXT). The SDK installs both under `share/alps/licenses/aleav2/`. Publication pointers in imported comments have been removed; the source relationships below refer to the remaining code.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `alea/include/alps/alea.hpp` | `include/alps/alea.hpp` | Copied |
| `alea/include/alps/alea/**` | `include/alps/alea/**` | Copied with the local changes listed below; `mpi.hpp` is installed only by MPI builds |
| `alea/src/*.cpp` (13 files) | `src/` | Copied except for the zero-step correction in `propagation.cpp` |
| `utilities/include/alps/common/view.hpp` | `include/alps/common/view.hpp` | Copied |
| `utilities/include/alps/serialization/{core,eigen,util}.hpp` | `include/alps/serialization/` | Copied, except that `core.hpp` gains an include (see below) |
| `alea/test/*.cpp` (10 files), `alea/test/dataset.hpp` | `tests/` | Copied with local regression cases in `transform.cpp` and `compare_test.cpp` |
| `utilities/include/alps/testing/near.hpp` | `tests/support/alps/testing/near.hpp` in the repository root | Adapted to enforce absolute matrix tolerances; shared test support |
| `COPYRIGHT.TXT`, `LICENSE.TXT` | here | Copied |

ALPSCore keeps `view.hpp` and the serialization headers in its utilities
module. In ALPS only ALEA v2 uses them, so they live here, at ALPSCore's include
spellings.

The CMake files, `tests/coexistence.cpp` and this README are ALPS's own.

### Local changes

`transformer.hpp` corrects an inherited dimension reversal: a rectangular matrix
maps its column count to its row count. The added `jacobian.rectangular` test
checks both expanding and contracting transformations, their Jacobians and
transformed means.

`propagation.cpp` supplies a positive, coordinate-scaled finite-difference step
when `jacobian` receives zero, as automatic propagation does for constant
observables. Nonzero steps retain their existing behavior. `propagation.hpp`
documents this extension. Regression tests check the Jacobian at zero and large
coordinates and preserve zero variance through identity transformations for
variance, covariance, autocorrelation and batch results.

The shared `near.hpp` helper now enforces an absolute error bound for each matrix
entry instead of a relative matrix norm. `compare_test.cpp` adds regression
coverage; the [test-support README](../../../tests/support/README.md) records
the adaptation.

Three imported headers relied on their includers for some declarations. Each
now includes what it needs, so it compiles on its own, and a comment marks the
addition:

- `alps/alea/convert.hpp` includes `mean.hpp`, `variance.hpp`, `covariance.hpp`
  and `batch.hpp`. They define the results it joins and give `var_result` and
  `cov_result` their default strategy arguments.
- `alps/alea/internal/outer.hpp` includes `<cassert>`; GCC rejects it otherwise.
- `alps/serialization/core.hpp` includes `<cstdint>`. On Linux, nothing else it
  includes declares `uint32_t` and `uint64_t`.

The test registration also differs:

- ALPS runs each GoogleTest case as its own CTest test, so the `twogauss` cases,
  which all write `twogauss.hdf5`, share a CTest resource lock.
- `mpi_twogauss` keeps its own `main()` and runs on 2 and 3 ranks. ALPSCore ran
  it on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- The tests' `TYPED_TEST_CASE` deprecation warnings are silenced.
