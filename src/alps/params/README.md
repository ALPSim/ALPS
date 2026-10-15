# Parameters

`ALPS::params` holds two parameter classes.

- ALPSCore's parameters, `alps::params_ns::params` and the `alps::params_ns::dictionary`
  it derives from, also named `alps::params` and `alps::dictionary`, in
  `<alps/params.hpp>`, `<alps/dictionary.hpp>` and `<alps/params/*.hpp>`. They
  are copied from the [ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore)
  at commit `7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25). They read
  INI files and command lines with the iniparser C library in `iniparser/`.
- The NGS parameters, `alps::params` in `<alps/ngs/params.hpp>`, with
  `<alps/ngs/detail/*.hpp>` and `src/ngs/`, are ALPS's own. The text/XML and
  older `Parameters` conversion adapters in `adapters/` belong to them and are
  compiled into `ALPS::alps`.

Both define `alps::params`, so a translation unit can include
`<alps/params.hpp>` or `<alps/ngs/params.hpp>`, not both. A program can use
both classes from separate translation units.

## Building and linking

```cmake
find_package(ALPS CONFIG REQUIRED COMPONENTS params)
target_link_libraries(app PRIVATE ALPS::params)
```

`ALPS::params` links `ALPS::utilities` and `ALPS::hdf5`, and MPI in MPI builds.
ALPSCore's parameters are stored through ALPS's `alps::hdf5::archive`, and in
MPI builds they are broadcast over ALPSCore's `<alps/utilities/mpi.hpp>`
communicator. iniparser is compiled into `ALPS::params`; its headers are not
installed.

## Provenance

The ALPSCore files are under the MIT license in [LICENSE.TXT](LICENSE.TXT).
Each refers to [COPYRIGHT.TXT](COPYRIGHT.TXT), [LICENSE.TXT](LICENSE.TXT) and
[ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT), which are ALPSCore's top-level files, copied
unchanged. The SDK installs them under `share/alps/licenses/params/`.

iniparser is by Nicolas Devillard, under the MIT license in
[iniparser/LICENSE](iniparser/LICENSE), and was modified for ALPSCore by
Alexander Gaenko, as `iniparser/src/iniparser.c` notes. The SDK installs its
`LICENSE` and `AUTHORS` under `share/alps/licenses/params/iniparser/`. Keep the
notices when modifying these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `params/include/alps/params.hpp`, `params/include/alps/dictionary.hpp` | `include/alps/` | Copied unchanged |
| `params/include/alps/params/*.hpp` (10 files) | `include/alps/params/` | Copied unchanged, except that three headers gain includes (see below) |
| `params/src/*.cpp` (4 files) | `src/` | Copied unchanged |
| `params/test/*.cpp` (22 files), `params_test_support.hpp`, `dict_values_test.hpp` | `tests/` | Copied unchanged |
| `common/deps/iniparser/` | `iniparser/` | Copied unchanged, without its `CMakeLists.txt` |
| `COPYRIGHT.TXT`, `LICENSE.TXT`, `ACKNOWLEDGE.TXT` | here | Copied unchanged |

Not copied: the retired tests in `params/test/old/`, the two scripts that
generated test code (`params_def_gen_test_helper.sh` and
`params_gen_equality_test_helper.py`), and ALPSCore's build files.

The utilities ALPSCore's parameters use are in `ALPS::utilities`; see the
[utilities README](../utilities/README.md). ALPSCore's `alps::testing::unique_file`,
which the tests use, is in [test support](../../../tests/support/README.md). The
CMake files and this README are ALPS's own.

### Local changes

Three imported headers relied on their includers for some declarations. Each
now includes what it needs, so it compiles on its own, and a comment marks the
addition:

- `alps/params/dict_exceptions.hpp` includes `<stdexcept>` and `<string>`.
- `alps/params/hdf5_variant.hpp` and `alps/params/mpi_variant.hpp` include
  `alps/params/dict_types.hpp`, which declares the value types they serialize.

`alps/params/dict_value_impl.hpp` and `alps/params/params_impl.hpp` are the
second halves of `dict_value.hpp` and `params.hpp`, which include them at
their end; they are not meant to be included on their own. The SDK header
contracts check the other headers.

### Test registration

- ALPS runs each GoogleTest case as its own CTest test, labelled `params`. A
  test whose name does not start with `params` gets the prefix: ALPSCore's
  `dictionary` is `params_dictionary`.
- The MPI tests (`serialize_variant_mpi`, `dictionary_bcast_mpi` and
  `params_bcast_mpi`) keep their own `main()` and run on 2 and 3 ranks.
  ALPSCore ran them on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- The tests' `TYPED_TEST_CASE` deprecation warnings are silenced.

## Citing

ALPSCore's [ACKNOWLEDGE.TXT](ACKNOWLEDGE.TXT) asks publications that use the ALPS
Core Libraries to cite M. Wallerberger et al., "Updated Core Libraries of the
ALPS Project", [arXiv:1811.08331](https://arxiv.org/abs/1811.08331) (2018),
alongside the ALPS paper. The ALPS citation catalog lists it under the `params`
component.
