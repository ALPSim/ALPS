# Shared test support

`alps_add_gtest` and `alps_add_mpi_gtest` put this directory on every test's
include path.

| File | Purpose |
| --- | --- |
| `mpi_main.cpp` | `main()` for MPI tests: starts MPI and combines the ranks' results |
| `alps/testing/temporary_directory.hpp` | A private working directory per test case |
| `alps/testing/stream_fixture.hpp` | Redirects `std::cin` and `std::cout` to compare a test's output with a reference file |

## Copied from ALPSCore

These are copied from ALPSCore commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969`. Their notices refer to ALPSCore's top-level
[`COPYRIGHT.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/COPYRIGHT.TXT), [`LICENSE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/LICENSE.TXT)
and [`ACKNOWLEDGE.TXT`](https://github.com/ALPSCore/ALPSCore/blob/7146b9e1f017938a94e5dae35d88467cc5ba7969/ACKNOWLEDGE.TXT).
They remain unchanged except for `near.hpp`: its matrix comparison now checks
the absolute error of every entry, matching its documented tolerance and scalar
overload, instead of using Eigen's relative `isApprox`. Regression coverage in
ALEA v2's `compare_test.cpp` includes large values, comparisons to zero, complex
entries, empty/mismatched shapes and nonfinite errors.

| ALPSCore path | Here |
| --- | --- |
| `utilities/include/alps/utilities/gtest_par_xml_output.hpp` | `alps/utilities/gtest_par_xml_output.hpp` |
| `utilities/include/alps/utilities/fs/{get_basename,get_dirname,get_extension,remove_extensions}.hpp` | `alps/utilities/fs/` |
| `utilities/src/{gtest_par_xml_output,filename_operations}.cpp` | `src/` |
| `utilities/include/alps/testing/near.hpp` | `alps/testing/near.hpp`, used by the ALEA tests |

`alps_add_mpi_gtest(... OWN_MAIN)` links `gtest_par_xml_output`, which ALPSCore's
MPI tests call from their `main()`. `alps/config.hpp` is ALPS's own: it stands in
for ALPSCore's generated configuration header by including `<alps/config.h>`.
