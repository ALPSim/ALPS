# Parameters

`ALPS::params` holds ALPSCore's parameters, `alps::params_ns::params` and the
`alps::params_ns::dictionary` it derives from, also named `alps::params` and
`alps::dictionary`, in `<alps/params.hpp>`, `<alps/dictionary.hpp>` and
`<alps/params/*.hpp>`. They are copied from the
[ALPS Core Libraries](https://github.com/ALPSCore/ALPSCore) at commit
`7146b9e1f017938a94e5dae35d88467cc5ba7969` (2026-09-25), and read INI files and
command lines with the iniparser C library in `iniparser/`. ALPS's programs,
solvers and pyalps use them. They replace the NGS parameters of
`<alps/ngs/params.hpp>`, which are removed.

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

## Using them

A program defines each parameter it reads, with its type, a description and,
where one applies, a default. It then checks that every parameter without a
default was given and parses:

```cpp
alps::params parameters = alps::make_parameters_from_file("param.txt");
parameters.define<int>("L", "number of sites")
          .define<double>("T", 1.0, "temperature");
if (parameters.has_missing(std::cerr))
    return EXIT_FAILURE;
const int length = parameters["L"].as<int>();
```

Reading a parameter that is not defined throws. `has_unused(std::ostream&)`
lists the parameters that were supplied but never defined; the solvers print
them as a warning. A Monte Carlo simulation puts its definitions in a static
`define_parameters()` that calls its base's: `alps::mcbase::define_parameters()`
defines `SEED` (default 42), and `alps::mcmpiadapter` adds `Tmin` and `Tmax`
(defaults 1 and 600), the bounds in seconds between its checks of the completed
fraction. The solvers and pyalps's simulation bindings call it before
constructing a simulation. In MPI programs, `alps::broadcast(communicator,
parameters)` from `<alps/mc/api.hpp>` takes an `alps::mpi::communicator`.

## Where values come from

Each source supplies text, as an INI file does, and `define<T>()` parses it with
the declared type. A value that does not parse is an error naming the
parameter, reported by `has_missing()`; so is a non-integral number for an
integer parameter.

| Source | Read with |
| --- | --- |
| INI file, command line | ALPSCore's constructors, `alps::params(argc, argv)` and `alps::params(filename)` |
| ALPS parameter file (`L=10; T=2.25;`) | `alps::make_parameters_from_file(path)`, `<alps/ngs/make_parameters_from_file.hpp>` |
| XML job file (`<PARAMETERS>`) | `alps::make_parameters_from_xml(path)`, `<alps/ngs/make_parameters_from_xml.hpp>`; `SEED` is `0` unless the file sets it |
| HDF5 group | `archive["/parameters"] >> parameters`; see the HDF5 local change below |
| Python dict | `pyalps.ngs.params(dict)` and the solver bindings; see [Python](#python) |

`alps::make_deprecated_parameters(parameters)`, in
`<alps/ngs/make_deprecated_parameters.hpp>`, converts to the legacy
`alps::Parameters` that the lattice and model libraries read. It passes the
supplied text of every parameter, including those the program never defines,
and leaves defaults to the code that reads `alps::Parameters`. The NGS
parameters passed numbers with six significant digits; the text is now passed
as supplied.

These three adapters, in `adapters/`, are ALPS's own and are compiled into
`ALPS::alps`.

## Python

`pyalps.ngs.params` wraps the same class. A value is converted to text when it
is stored, so the params keep no reference to the Python object. Booleans
become `1` and `0`, integers keep every digit, and sequences become
comma-separated lists. `None`, complex numbers, multidimensional arrays, nested
sequences and dicts, and list elements containing a comma are a `TypeError`
naming the parameter.

A parameter that a program has defined reads back with its type. One that is
only supplied reads back from its text, as `supplied_value()` in
`<alps/params/text_values.hpp>` reads it: an integer, a float, a NumPy array of
numbers, or else the string. Until a program defines them, `"42"` reads back as
`42`, `2.0` as `2`, `[5]` as `5` and `[True, False]` as `array([1, 0])`. Keys
iterate in sorted order. `del parameters[key]` rebuilds the params from their
values, which drops every definition.

## Provenance

The ALPSCore files are under the MIT license in [LICENSE.TXT](LICENSE.TXT).
The imported copyright and license notices are retained in [COPYRIGHT.TXT](COPYRIGHT.TXT) and [LICENSE.TXT](LICENSE.TXT). The SDK installs both under `share/alps/licenses/params/`. Publication pointers have been removed from imported comments.

iniparser is by Nicolas Devillard, under the MIT license in
[iniparser/LICENSE](iniparser/LICENSE), and was modified for ALPSCore by
Alexander Gaenko, as `iniparser/src/iniparser.c` notes. The SDK installs its
`LICENSE` and `AUTHORS` under `share/alps/licenses/params/iniparser/`. Keep the
notices when modifying these files.

| ALPSCore path | Here | Relationship |
| --- | --- | --- |
| `params/include/alps/params.hpp`, `params/include/alps/dictionary.hpp` | `include/alps/` | Copied, except that `params.hpp` gains two members (see below) |
| `params/include/alps/params/*.hpp` (10 files) | `include/alps/params/` | Copied, except that three headers gain includes (see below) |
| `params/src/*.cpp` (4 files) | `src/` | Copied, except `params.cpp` (see below) |
| `params/test/*.cpp` (22 files), `params_test_support.hpp`, `dict_values_test.hpp` | `tests/` | Copied |
| `common/deps/iniparser/src/`, `LICENSE`, `AUTHORS` | `iniparser/` | Copied; upstream build instructions and FAQs for files outside this integration are omitted |
| `COPYRIGHT.TXT`, `LICENSE.TXT` | here | Copied |

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

Two additions let other sources supply values the way INI files and the
command line do, and let programs read what was supplied. Comments marked
`ALPS:` in `params.hpp` and `params.cpp` mark them, and `tests/params_supply.cpp`
tests them:

- `params::supply(name, text)` stores text for a parameter that is not yet
  defined. `define<T>()` parses it as it parses a value from a file.
- `params::supplied_text()` returns the supplied text of every supplied
  parameter, by name, including parameters the program never defines.

When an integer parameter rejects whole-number text such as `24.0` or `1e6`,
the error in `params_impl.hpp` says how to write it: `write 1000000`. A
comment marked `ALPS:` marks the addition.

`params::load()` also reads HDF5 groups without ALPSCore's attributes, as
pyalps, `alps::Parameters::save` and the NGS parameters wrote `/parameters`. It
supplies each value as text, as an INI file would, so `define<T>()` parses it
with the type the program declares; a typed load would instead fail on the
first value stored with another type, such as a Python integer for a `double`.
ALPSCore's own groups, recognized by their `@status` attribute, load as
before. The text comes from `<alps/params/text_values.hpp>`
(`src/text_values.cpp`), which is ALPS's own, and `tests/params_plain_group.cpp`
tests it:

| Stored value | Supplied text |
| --- | --- |
| Integer of any width | Decimal: `1000` |
| Boolean, stored as `int8` | `1` or `0`, which parse as `bool` and as integers |
| Floating-point number | Integral values without a fraction: `24`, `1000000`; others with the digits they need to read back unchanged: `0.123456789` |
| String | Unchanged |
| One-dimensional array, or a group of scalars `0`, `1`, ... (a Python list) | Elements separated by commas: `4,4` |
| Complex number, multidimensional array or other group | An error naming the parameter |

An integer parameter therefore accepts `24.0` and `1e6` stored as numbers,
but not `2.7`, nor the strings `"4.0"` or `"1e6"`: these are errors naming the
parameter.

`params::save()` writes a parameter that is supplied but not defined as a
value too, so `/parameters` keeps every parameter as a dataset, as ALPS's
other writers left it, and pyalps's readers, which list the group's children,
find them. ALPSCore's `save()` keeps such a parameter only in the group's
`@ini_values` attribute. The value is `supplied_value()` of the text, from
`<alps/params/text_values.hpp>`: an integer that fits in a `long`, another
number as a `double`, a list of numbers as a vector, or else the string. Each
carries an `@supplied_only` attribute, and `load()` drops these values, so
the parameter stays supplied text that `define<T>()` parses with the
program's type. `save()` also creates the group when there is nothing to
write, where ALPSCore's failed writing its attributes.

### Test registration

- ALPS runs each GoogleTest case as its own CTest test, labelled `params`. A
  test whose name does not start with `params` gets the prefix: ALPSCore's
  `dictionary` is `params_dictionary`.
- The MPI tests (`serialize_variant_mpi`, `dictionary_bcast_mpi` and
  `params_bcast_mpi`) keep their own `main()` and run on 2 and 3 ranks.
  ALPSCore ran them on `ALPS_TEST_MPI_NPROC` ranks, 1 by default.
- The tests' `TYPED_TEST_CASE` deprecation warnings are silenced.
