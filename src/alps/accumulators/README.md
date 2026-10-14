# Accumulators

This directory holds two accumulator libraries.

- `ALPS::accumulators` is the accumulators library of the
  [ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore): `alps::accumulators`
  with `<alps/accumulators.hpp>` and `<alps/accumulators/*.hpp>`. It is copied
  from ALPSCore at commit `7146b9e1f017938a94e5dae35d88467cc5ba7969`
  (2026-09-25).
- The NGS accumulators, `alps::accumulator` with
  `<alps/ngs/accumulator/*.hpp>`, are ALPS's own. They are compiled into
  `ALPS::alps` only with `ALPS_NGS_USE_NEW_ALEA=ON`, and their tests in
  `tests/legacy_alea/` and `tests/legacy_accumulator/` are not built.

The two use different namespaces and include paths, with one overlap: both
read the macro `ALPS_ACCUMULATOR_VALUE_TYPES`. `<alps/config.hpp>` defines it as
ALPSCore's list, and the NGS headers define their own list only if it is not
yet defined. A translation unit that includes both therefore gives the NGS
accumulators ALPSCore's list, or, if the NGS headers come first, a macro
redefinition warning.

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
Every imported file refers to [COPYRIGHT.TXT](COPYRIGHT.TXT),
[LICENSE.TXT](LICENSE.TXT) and [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT), which are
ALPSCore's top-level files, copied unchanged. The SDK installs them under
`share/alps/licenses/accumulators/`. Keep the notices when modifying these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `accumulators/include/alps/accumulators.hpp` | `include/alps/accumulators.hpp` | Copied unchanged |
| `accumulators/include/alps/accumulators/**` | `include/alps/accumulators/**` | Copied unchanged, except that `convergence.hpp` and `feature/binning_analysis.hpp` gain includes (see below) |
| `accumulators/src/*.cpp`, `accumulators/src/feature/*.cpp` (11 files) | `src/`, `src/feature/` | Copied unchanged |
| `accumulators/test/*.cpp` (31 files), `accumulators/test/accumulator_generator.hpp` | `tests/` | Copied unchanged, except two include lines (see below) |
| `utilities/test/vector_comparison_predicates.hpp` | `tests/vector_comparison_predicates.hpp` | Copied unchanged |
| `utilities/include/alps/testing/unique_file.hpp`, `utilities/src/unique_file.cpp` | `tests/support/` in the repository root | Copied unchanged; shared test support |
| `COPYRIGHT.TXT`, `LICENSE.TXT`, `ACKNOWLEDGE.TXT` | here | Copied unchanged |

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

## Citing

ALPSCore's [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT) asks publications that use the ALPS
Core Libraries to cite M. Wallerberger et al., "Updated Core Libraries of the
ALPS Project", [arXiv:1811.08331](https://arxiv.org/abs/1811.08331) (2018),
alongside the ALPS paper. The ALPS citation catalog lists it under the
`accumulators` component.
