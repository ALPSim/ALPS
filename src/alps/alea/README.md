# ALEA

`ALPS::alea` is the ALEA library of the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore): accumulators for
mean, variance, covariance, autocorrelation and batch statistics, error
propagation, and result serialization. It is copied from ALPSCore at commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25). Its include spellings
are ALPSCore's: `<alps/alea.hpp>` and `<alps/alea/*.hpp>`, in namespace
`alps::alea`.

The older ALEA in [legacy_alea/](../legacy_alea/) keeps `<alps/alea.h>` and
`<alps/alea/*.h>`. The two share the include prefix and namespace without
sharing a file name, and one program can use both; `alea_legacy_coexistence`
checks this, including both kinds of result in one archive.

## Building and linking

ALEA needs [Eigen](https://eigen.tuxfamily.org) 3.3 or newer; it was tested with
Eigen 5.0.1. ALPS builds `ALPS::alea` only when CMake finds Eigen, and reports
otherwise during configuration. Set `CMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON` to
leave it out. An installed SDK provides it as a separate component, so other
consumers do not need Eigen:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS alea)
target_link_libraries(app PRIVATE ALPS::alea)
```

ALEA stores results through ALPS's `alps::hdf5::archive`
([`hdf5_serializer`](include/alps/alea/hdf5.hpp)).

## Provenance

ALPSCore was licensed under the GPL until commit `401fdd50` (2025-09-14,
"Replace GPL with permissive license terms"). The imported revision is under
the MIT license in [LICENSE.TXT](LICENSE.TXT). Every imported file refers to
[COPYRIGHT.TXT](COPYRIGHT.TXT), [LICENSE.TXT](LICENSE.TXT) and
[ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT), which are ALPSCore's top-level files, copied
unchanged. The SDK installs them under `share/alps/licenses/alea/`. Keep the
notices when modifying these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `alea/include/alps/alea.hpp` | `include/alps/alea.hpp` | Copied; see local changes |
| `alea/include/alps/alea/**` | `include/alps/alea/**` | Copied unchanged, except `mpi.hpp` |
| `alea/src/*.cpp` (13 files) | `src/` | Copied unchanged |
| `utilities/include/alps/common/view.hpp` | `include/alps/common/view.hpp` | Copied unchanged |
| `utilities/include/alps/serialization/{core,eigen,util}.hpp` | `include/alps/serialization/` | Copied unchanged |
| `alea/test/*.cpp` (9 files), `alea/test/dataset.hpp` | `tests/` | Copied unchanged, except `mpi_twogauss.cpp` |
| `utilities/include/alps/testing/near.hpp` | `tests/support/alps/testing/near.hpp` | Copied unchanged; used only by the tests |
| `COPYRIGHT.TXT`, `LICENSE.TXT`, `ACKNOWLEDGE.TXT` | here | Copied unchanged |

ALPSCore keeps `view.hpp` and the serialization headers in its utilities
module. In ALPS only ALEA uses them, so they live here, at ALPSCore's include
spellings.

The CMake files, `tests/legacy_coexistence.cpp` and this README are ALPS's own.

### Local changes

- `alps/alea.hpp` no longer includes `alps/alea/mpi.hpp`, which is not imported:
  it needs ALPSCore's `alps/utilities/mpi.hpp`, which ALPS does not have. For
  the same reason the MPI test `mpi_twogauss` is not ported. Without this change, an MPI-enabled ALPS build, which defines
  `ALPS_HAVE_MPI`, could not include the umbrella header.
- The test registration differs, not the test files. ALPS runs each GoogleTest
  case as its own CTest test, so the `twogauss` cases, which all write
  `twogauss.hdf5`, share a CTest resource lock. The tests' `TYPED_TEST_CASE`
  deprecation warnings are silenced.

## Citing

ALPSCore's [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT) asks publications that use the ALPS
Core Libraries to cite M. Wallerberger et al., "Updated Core Libraries of the
ALPS Project", [arXiv:1811.08331](https://arxiv.org/abs/1811.08331) (2018),
alongside the ALPS paper. The ALPS citation catalog lists it under the `alea`
component.
