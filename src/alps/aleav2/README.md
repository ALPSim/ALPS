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
otherwise during configuration. Set `CMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON` to
leave it out. An installed SDK provides it as a separate component, so other
consumers do not need Eigen:

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
Every imported file refers to
[COPYRIGHT.TXT](COPYRIGHT.TXT), [LICENSE.TXT](LICENSE.TXT) and
[ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT), which are ALPSCore's top-level files, copied
unchanged. The SDK installs them under `share/alps/licenses/aleav2/`. Keep the
notices when modifying these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `alea/include/alps/alea.hpp` | `include/alps/alea.hpp` | Copied unchanged |
| `alea/include/alps/alea/**` | `include/alps/alea/**` | Copied unchanged; `mpi.hpp` is installed only by MPI builds |
| `alea/src/*.cpp` (13 files) | `src/` | Copied unchanged |
| `utilities/include/alps/common/view.hpp` | `include/alps/common/view.hpp` | Copied unchanged |
| `utilities/include/alps/serialization/{core,eigen,util}.hpp` | `include/alps/serialization/` | Copied unchanged |
| `alea/test/*.cpp` (10 files), `alea/test/dataset.hpp` | `tests/` | Copied unchanged |
| `utilities/include/alps/testing/near.hpp` | `tests/support/alps/testing/near.hpp` in the repository root | Copied unchanged; shared test support |
| `COPYRIGHT.TXT`, `LICENSE.TXT`, `ACKNOWLEDGE.TXT` | here | Copied unchanged |

ALPSCore keeps `view.hpp` and the serialization headers in its utilities
module. In ALPS only ALEA v2 uses them, so they live here, at ALPSCore's include
spellings.

The CMake files, `tests/coexistence.cpp` and this README are ALPS's own.

### Local changes

No imported file is changed. The test registration differs:

- ALPS runs each GoogleTest case as its own CTest test, so the `twogauss` cases,
  which all write `twogauss.hdf5`, share a CTest resource lock.
- `mpi_twogauss` keeps its own `main()` and runs on 2 and 3 ranks. ALPSCore ran
  it on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- The tests' `TYPED_TEST_CASE` deprecation warnings are silenced.

## Citing

ALPSCore's [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT) asks publications that use the ALPS
Core Libraries to cite M. Wallerberger et al., "Updated Core Libraries of the
ALPS Project", [arXiv:1811.08331](https://arxiv.org/abs/1811.08331) (2018),
alongside the ALPS paper. The ALPS citation catalog lists it under the `aleav2`
component.
