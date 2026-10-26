# Utilities

`ALPS::utilities` holds the general helpers every other module may use. Its
headers sit under `alps/utilities/`, the spelling ALPSCore uses, plus the
topical prefix `alps/type_traits/`. The directory stays flat.

A header belongs in `alps/utilities/` only if

1. more than one module uses it; otherwise it lives with its user,
2. nothing is planned to replace it; otherwise it stays with the legacy code, and
3. Boost idioms with a standard equivalent were replaced when it moved, without
   changing public signatures. Signature changes, such as
   `boost::filesystem::path` to `std::filesystem::path`, are later work.

Status *current* means the header needs no planned changes. *Carried over*
means it was moved as it was, and the note says what modernization remains.

| Header | Purpose | Status |
| --- | --- | --- |
| `cast.hpp` | `alps::cast<U>(T)` conversions, including strings and paths | current; ALPSCore name. ALPS keeps its filesystem hook and range-checked integer parsing |
| `stacktrace.hpp` | `alps::stacktrace()` and `ALPS_STACKTRACE` | current; ALPSCore name |
| `signal.hpp` | `alps::signal`: POSIX signal capture and the segfault handler | current; ALPSCore name. `set_fatal_cleanup()` is an ALPS extension used by `ALPS::hdf5` |
| `stringify.hpp` | `ALPS_STRINGIFY` | current; ALPSCore name |
| `short_print.hpp` | Abbreviated printing of values and sequences | current; ALPSCore name. ALPS adds `alps::multi_array` |
| `remove_cvr.hpp` | `alps::detail::remove_cvr` | carried over; ALPSCore name, superseded by `std::remove_cvref_t` in C++20 |
| `type_wrapper.hpp` | `alps::detail::type_wrapper` | current; ALPSCore name |
| `size.hpp` | `alps::size()` for scalars and sequences | current; ALPSCore name |
| `temporary_filename.hpp` | Private temporary file names | current; ALPSCore name. ALPS keeps its by-value signature |
| `make_copy.hpp` | `alps::make_copy()` | current; ALPSCore name |
| `mpi.hpp` | `alps::mpi`: ALPSCore's MPI communicator, environment, broadcast and reductions | current; copied from ALPSCore with local adaptations, MPI builds only |
| `boost_mpi.hpp` | `alps::mpi::reduce` for `std::vector`, on ALPSCore's communicator; despite the name, no Boost.MPI | current; copied from ALPSCore, MPI builds only |
| `mpi_vector.hpp`, `mpi_pair.hpp`, `mpi_map.hpp` | `alps::mpi::broadcast` for `std::vector`, `std::pair` and `std::map` | current; copied from ALPSCore, MPI builds only |
| `fs/get_basename.hpp`, `fs/get_dirname.hpp`, `fs/get_extension.hpp`, `fs/remove_extensions.hpp` | File-name operations on strings | current; copied from ALPSCore with `src/filename_operations.cpp` |
| `deprecated.hpp` | `ALPS_DEPRECATED` | current; copied from ALPSCore, and defines the macro only if ALPS's `<alps/export.h>` has not |
| `copyright.hpp` | Program banners, version and license text | current |
| `citations.hpp` | Citation notices from the generated catalog | current |
| `os.hpp` | Host, user and directory queries | carried over; returns `boost::filesystem::path` |
| `vmusage.hpp` | Process memory usage | carried over; Linux `/proc` only |
| `vectorio.hpp` | Reading and writing whitespace-separated vectors | carried over; Boost stream idioms |
| `encode.hpp` | HDF5 name escaping | carried over; Boost Software License notice, `boost::lexical_cast` |
| `resize.hpp` | `alps::resize_same_as()` for containers and multi-arrays | carried over; used by `numerics` and ALEA |
| `data.hpp` | `alps::data()` for contiguous containers | carried over; overlaps `std::data` |

`alps/testing/fp_compare.hpp`, ALPSCore's floating-point comparison
`alps::testing::is_near`, is installed with these headers as in ALPSCore, where
library code (`ALPS::params`) uses it too. Its publication pointer has been removed.

Not yet moved: `alps/ngs/{config,hash,mutex,sleep,thread_exceptions}.hpp`, the
top-level `alps/cctype.h`, `alps/factory.h`, `alps/lambda.hpp` and
`alps/stringvalue.h`, and the type traits ALPSCore lacks. `alps/ngs/boost_mpi.hpp`
stays where it is: ALPSCore's `alps/utilities/boost_mpi.hpp` is a different API.

[`ALPS::gf`](../gf/README.md) installs four of ALPSCore's `alps/type_traits/`
headers. Only GF uses them, so by rule 1 they live with it.

Old spellings of the moved headers forward from the [legacy layer](../../../legacy/README.md)
while `ALPS_LEGACY_API` is on.

## ALPSCore's MPI wrapper

`mpi.hpp` is copied from ALPSCore commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969`
(`utilities/include/alps/utilities/mpi.hpp`) with the macro adaptation marked by
comments: ALPSCore's unprefixed macro `EXCEPTION_CATCH` is renamed
`ALPS_MPI_EXCEPTION_CATCH` and undefined at the end of the header, so it no
longer leaks into every file that includes it, directly or through ALEA v2.
Its notice refers to ALPSCore's top-level
[`COPYRIGHT.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/COPYRIGHT.TXT),
[`LICENSE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/LICENSE.TXT).
Publication pointers in imported comments have been removed.
It includes `<mpi.h>`, so only MPI builds install it. In ALPSCore, ALEA, the
Green's functions, parameters and Monte Carlo code all use it, so it lives here
rather than with one module. `boost_mpi.hpp`, from the same commit and copied
with publication pointers removed, adds its `std::vector` reductions, which ALPSCore's accumulators use.

The container broadcasts `mpi_vector.hpp`, `mpi_pair.hpp` and `mpi_map.hpp`,
which ALPSCore's parameters use, are copied from the same commit.

Its `environment`, like Boost.MPI's, leaves MPI alone if something else
initialized it.

`tests/` holds ALPSCore's tests of the header (`utilities/test/`):
`mpi_utils_env`, `mpi_utils_comm`, `mpi_utils_bcast`, `mpi_utils_reduce` and the
two `mpi_utils_env_exceptions`, with `mpi_test_support.{cpp,hpp}` and
`test_utils.hpp`, with publication pointers removed from their comments. They run on 2 and 3 ranks; ALPSCore ran them on
`ALPS_TEST_MPI_NPROC` ranks, 1 by default. `mpi_utils_bcast_optional` is not
ported, since it tests `mpi_optional.hpp`, which ALPS does not have.
