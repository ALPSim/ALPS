# GF

`ALPS::gf` is the Green's function library of the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore). It provides Green's
function containers on Matsubara, imaginary-time, Legendre, Chebyshev, numerical
and index meshes, with high-frequency tails, Fourier transforms, piecewise
polynomials, HDF5 storage and MPI broadcast. It is copied from ALPSCore at
commit `7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25). Its include
spellings are ALPSCore's: `<alps/gf/gf.hpp>` and `<alps/gf/*.hpp>`, in namespace
`alps::gf`. The tensors it stores data in come with it: `<alps/numeric/tensors.hpp>`,
in namespace `alps::numerics`.

## Building and linking

GF needs [Eigen](https://eigen.tuxfamily.org) 3.3 or newer, like
[ALEA v2](../aleav2/README.md); it was tested with Eigen 5.0.1. ALPS builds `ALPS::gf`
only when CMake finds Eigen, and reports otherwise during configuration. Set
`CMAKE_DISABLE_FIND_PACKAGE_Eigen3=ON` to leave it out. An installed SDK
provides it as a separate component:

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS gf)
target_link_libraries(app PRIVATE ALPS::gf)
```

Include `<alps/gf/gf.hpp>`, or a header beside it such as
`<alps/gf/mesh.hpp>`. As in ALPSCore, `<alps/gf/mesh/*.hpp>` and
`<alps/gf/mpi_bcast.hpp>` are parts of those and do not compile alone.

GF stores Green's functions through ALPS's `alps::hdf5::archive`. In an MPI
build, `ALPS::gf` links `ALPS::utilities` and MPI, and installs
`<alps/gf/mpi_bcast.hpp>`. gf's `broadcast()` takes the `alps::mpi::communicator`
of ALPSCore's [`<alps/utilities/mpi.hpp>`](../utilities/README.md#alpscores-mpi-wrapper),
which `ALPS::utilities` installs in MPI builds. GF includes both only in MPI
builds.

`<alps/gf/mesh.hpp>` warns unless `BOOST_DISABLE_ASSERTS` and `NDEBUG` are
defined, or `ALPS_GF_DEBUG` is. ALPSCore's build defined `BOOST_DISABLE_ASSERTS`
for optimized configurations and `ALPS_GF_DEBUG` for Debug, but only for its own
code. ALPS does the same for `ALPS::gf` and its tests. Without
`BOOST_DISABLE_ASSERTS`, a broadcast also checks that all ranks agree on the
array dimensions, which costs extra communication. A consumer can match
ALPSCore's build with:

```cmake
target_compile_definitions(app PRIVATE
  "$<IF:$<CONFIG:Debug>,ALPS_GF_DEBUG,BOOST_DISABLE_ASSERTS>")
```

With libc++, ALPS's own build warns that `gf_base.hpp` ignores the result of
`std::tie` (`-Wunused-result`) where Green's functions are assigned, saved or
loaded. The code uses `std::tie` to expand a parameter pack; the warning is
harmless, and the copy keeps ALPSCore's code. SDK consumers do not see it,
because CMake includes an installed SDK's headers as system headers.

## Provenance

The imported revision is under the MIT license in [LICENSE.TXT](LICENSE.TXT).
The imported copyright and license notices are retained in [COPYRIGHT.TXT](COPYRIGHT.TXT) and [LICENSE.TXT](LICENSE.TXT). The SDK installs both under `share/alps/licenses/gf/`. Publication pointers in imported comments have been removed; the source relationships below refer to the remaining code.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `gf/include/alps/gf/**` (18 headers) | `include/alps/gf/**` | Copied; `mpi_bcast.hpp` is installed only by MPI builds |
| `gf/src/mesh.cpp` | `src/mesh.cpp` | Copied |
| `gf/test/*.cpp`, `gf/test/*.hpp` (30 files) | `tests/` | Copied |
| `utilities/include/alps/numeric/tensors.hpp`, `utilities/include/alps/numeric/tensors/*.hpp` (4 headers) | `include/alps/numeric/` | Copied |
| `utilities/include/alps/type_traits/{are_all_integrals,common,index_sequence,tuple_traits}.hpp` | `include/alps/type_traits/` | Copied, except that `index_sequence.hpp` and `tuple_traits.hpp` gain an include (see below) |
| `hdf5/include/alps/hdf5/tensor.hpp` | `include/alps/hdf5/tensor.hpp` | Copied |
| `utilities/include/alps/testing/unique_file.hpp` | `tests/support/alps/testing/unique_file.hpp` | Copied |
| `utilities/src/unique_file.cpp` | `tests/support/src/unique_file.cpp` | Copied |
| `COPYRIGHT.TXT`, `LICENSE.TXT` | here | Copied |

ALPSCore keeps the tensors and type traits in its utilities module and
`tensor.hpp` in its HDF5 module. In ALPS only GF uses them, so they live here, at
ALPSCore's include spellings. `tensor.hpp` throws ALPSCore's
`alps::hdf5::wrong_dimensions`, which ALPS's `<alps/hdf5/errors.hpp>` now
defines too.

The tests also use ALPSCore's `near.hpp`, its MPI test helpers and the
`<alps/config.hpp>` stand-in from the repository's
[shared test support](../../../tests/support/README.md). The CMake files,
This README is ALPS's own.

### Local changes

Two imported type traits relied on their includers for some declarations. Each
now includes what it needs, so it compiles on its own, and a comment marks the
addition:

- `alps/type_traits/index_sequence.hpp` includes `<cstddef>` for `size_t`.
- `alps/type_traits/tuple_traits.hpp` includes `index_sequence.hpp`, whose
  `alps::index_sequence` it uses.

The test registration also differs:

- ALPS runs each GoogleTest case as its own CTest test. The cases that write a
  fixed file name each use their own, so they run in parallel.
- The MPI tests keep their own `main()` and run on 2 and 3 ranks. ALPSCore ran
  them on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- `gf_check_eigen` has no test cases, so it is compiled but not run.
- `gf1_print_test`'s `TYPED_TEST_CASE` deprecation warning is silenced.
