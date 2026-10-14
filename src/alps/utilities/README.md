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
| `mpi.hpp` | `alps::mpi`: ALPSCore's MPI communicator, environment, broadcast and reductions | current; copied unchanged from ALPSCore, MPI builds only |
| `copyright.hpp` | Program banners, version and license text | current |
| `citations.hpp` | Citation notices from the generated catalog | current |
| `os.hpp` | Host, user and directory queries | carried over; returns `boost::filesystem::path` |
| `vmusage.hpp` | Process memory usage | carried over; Linux `/proc` only |
| `vectorio.hpp` | Reading and writing whitespace-separated vectors | carried over; Boost stream idioms |
| `encode.hpp` | HDF5 name escaping | carried over; Boost Software License notice, `boost::lexical_cast` |
| `resize.hpp` | `alps::resize_same_as()` for containers and multi-arrays | carried over; used by `numerics` and ALEA |
| `data.hpp` | `alps::data()` for contiguous containers | carried over; overlaps `std::data` |

Not yet moved: `alps/ngs/{config,hash,mutex,sleep,thread_exceptions}.hpp`, the
top-level `alps/cctype.h`, `alps/factory.h`, `alps/lambda.hpp` and
`alps/stringvalue.h`, and the type traits ALPSCore lacks. `alps/ngs/boost_mpi.hpp`
stays where it is: ALPSCore's header of the same name is a different API.

Old spellings of the moved headers forward from the [legacy layer](../../../legacy/README.md)
while `ALPS_LEGACY_API` is on.

## ALPSCore's MPI wrapper

`mpi.hpp` is copied unchanged from ALPSCore commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969`
(`utilities/include/alps/utilities/mpi.hpp`). Its notice refers to ALPSCore's
top-level
[`COPYRIGHT.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/COPYRIGHT.TXT),
[`LICENSE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/LICENSE.TXT)
and
[`ACKNOWLEDGE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/ACKNOWLEDGE.TXT).
It includes `<mpi.h>`, so only MPI builds install it. In ALPSCore, ALEA, the
Green's functions, parameters and Monte Carlo code all use it, so it lives here
rather than with one module.

Its names share namespace `alps::mpi` with ALPS's `<alps/ngs/mpi.hpp>`, without
overlapping; `utilities_mpi_coexistence` uses both in one program. Its
`environment`, like Boost.MPI's, leaves MPI alone if something else initialized
it.

`tests/` holds ALPSCore's tests of the header (`utilities/test/`):
`mpi_utils_env`, `mpi_utils_comm`, `mpi_utils_bcast`, `mpi_utils_reduce` and the
two `mpi_utils_env_exceptions`, with `mpi_test_support.{cpp,hpp}` and
`test_utils.hpp`, all unchanged. They run on 2 and 3 ranks; ALPSCore ran them on
`ALPS_TEST_MPI_NPROC` ranks, 1 by default. `mpi_utils_bcast_optional` is not
ported, since it tests `mpi_optional.hpp`, which ALPS does not have.
