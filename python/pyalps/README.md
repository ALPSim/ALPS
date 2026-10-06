# pyalps

Python applications and libraries for the Algorithms and Libraries for
Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install pyalps
```

Keep this environment separate from a source-installed ALPS SDK. Pip and CMake
do not share file ownership: installing a bundled wheel into the SDK's prefix
can overwrite its commands, and uninstalling the wheel will not restore them.
For Python bindings that share an SDK prefix, use the bindings-only build below.

## Command-line applications

The wheel installs shell commands for all bundled applications, including
`spinmc`, `loop`, `worm`, `dmrg`, and `sparsediag`, plus the tutorial tools
`parameter2xml`, `printgraph`, `convert2xml`, `convert2text`, `plot2text`,
`plot2gp`, `plot2xmgr`, `snap2vtk`, and `maxent`. Activate your Python environment so its
`bin` directory is on `PATH`, then follow the command-line tutorials:

```sh
parameter2xml simulation.in
spinmc simulation.in.in.xml
printgraph lattice.in
```

`parameter2xml simulation.in` creates `simulation.in.in.xml` (the job file)
and the task input files. Download tutorial inputs from the
[ALPS tutorial pages](https://alps.comp-phys.org/tutorials/); tutorial
directories are not included in the wheel.

The launchers execute the binaries bundled with this Python installation,
with the bundled XML library as the default `ALPS_XML_PATH`. An explicit
`ALPS_XML_PATH` override is respected. Pip does not edit shell startup files;
for a user installation, ensure Python's scripts directory is on `PATH`.
When an environment is active, its ALPS commands may take precedence over
another ALPS installation. Use a full executable path to select that install.

The four text/plot exporters use the packaged XSL stylesheets through `lxml`,
which pip installs automatically; no system `xsltproc` is needed. They take
an XML filename (or `-` for standard input) and write to standard output:

```sh
convert2text simulation.in.task1.out.xml > results.txt
plot2text susceptibility.plot.xml > susceptibility.txt
```

`convert2xml` converts legacy simulation checkpoints to XML, `snap2vtk`
converts simulation snapshots for visualization, and `maxent` runs analytic
continuation from an HDF5 parameter file. Gnuplot, Grace, and VTK viewers
are separate applications; the exporters produce their input files.

The obsolete `alpspython` wrapper has been removed: run `python` in your
active environment instead. `plot2mpl`, `extractmpl`, and their Python 2
code-generating stylesheet have also been removed; use `pyalps.plot` for
Matplotlib plots. The unused `transformall` helper and deprecated no-op
`changestylesheet` command have been removed. Other specialized SDK tools
remain available from source/SDK installations. `dirloop_sse_evaluate`
is not currently built and is not included in wheels.

Bundled applications are serial; installing `pyalps[mpi]`
does not turn them into MPI-enabled executables.

## Python and source builds

Install `pyalps[plot]` to use the Matplotlib plotting helpers.
Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

The bindings are built as a standalone `scikit-build-core` project using
nanobind. A source build requires Python 3.10 or newer, CMake 3.22 or newer,
Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed ALPS C++ SDK.
Point `ALPS_DIR` at the SDK's `share/alps` package directory.

The `wheel-deps` CMake preset builds the SDK exactly as the wheel CI does.
From the repository root:

```sh
cmake --preset wheel-deps
cmake --build --preset wheel-deps

ALPS_DIR="$PWD/_build/wheel-deps/install/share/alps" \
  python -m build --wheel python/pyalps
```

The wheel is written to `python/pyalps/dist` and can be installed
with `python -m pip install`. With ccache installed, configure with
`cmake --preset wheel-deps -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set
`CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to
speed up rebuilds.

`PYALPS_BUILD_APPLICATIONS=ON` is the default and preserves the MaxEnt,
CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.

`PYALPS_BUNDLE_APPLICATIONS=ON` is the default and copies the ALPS
application executables (`spinmc`, `dmrg`, `sparsediag`, `loop`, `qwl`, ...)
and the `parameter2xml`, `printgraph`, `convert2xml`, `snap2vtk`, and `maxent`
tools from the SDK into `pyalps/bin`,
together with the SDK's shared libraries in `pyalps/lib` that their `../lib` RPATH
resolves against. Both the shell launchers and `pyalps.runApplication('spinmc', ...)`
resolve that bundled executable directly; another command on `PATH` cannot
override it. The `wheel-deps` preset therefore builds the applications.

For a bindings-only installation, set the **environment variable**
`PYALPS_BUNDLE_APPLICATIONS=OFF` for the entire Python build:

```sh
ALPS_DIR="/path/to/sdk/share/alps" PYALPS_BUNDLE_APPLICATIONS=OFF \
  python -m pip install ./python/pyalps
```

This installs no command launchers, so the SDK retains ownership of its
executables even if it shares the Python environment's prefix. Add the SDK's
`bin` directory to `PATH` to use its shell commands. Python's application and
evaluation helpers use `ALPS_BIN_PATH` if explicitly set, otherwise the SDK
recorded at build time. A missing selected executable is an error; the helpers
do not silently select another installation from `PATH`.

The environment setting replaces the CMake-only
`-DPYALPS_BUNDLE_APPLICATIONS=OFF` option: wheel metadata is prepared before
CMake runs, so both stages must receive the same choice. Bindings-only wheels
remain tied to an installed SDK and are intended for local use.

In either build mode, pass a full executable path (or a relative path containing
a directory) to select a particular application, including a source-built MPI
application. For bare command names, bundled wheels always use their own
payload and ignore `ALPS_BIN_PATH`; bindings-only builds honor the SDK selection
described above. `ALPS_XML_PATH` remains an independent resource override.

## Free-threading and stable-ABI policy

pyalps ships per-version wheels (CPython 3.10–3.14) and deliberately opts
into neither of nanobind's special ABI modes:

- **Free-threading (3.13t/3.14t):** the extension modules do not declare
  free-threading support, so importing pyalps on a free-threaded
  interpreter re-enables the GIL for the process. That is intentional:
  the ALPS C++ library relies on the GIL as its lock around shared state
  (`mcobservable`'s reference-count table, the `alps::ngs::signal`
  singleton, `mcdata`'s lazily-computed statistics). Do not add
  `FREE_THREADED` to `nanobind_add_module` without first making that
  state thread-safe.
- **Stable ABI (abi3):** not enabled. Nanobind isolates stable-ABI and
  ordinary extensions from each other. ALPS supports downstream nanobind
  modules that derive from pyalps types, so an abi3 pyalps wheel would force
  every such consumer to use the limited API too. Per-version wheels preserve
  ordinary downstream extension interoperability.

## Compatibility and checkpoints

Parameters created from Python retain their Python values. NumPy arrays keep
array arithmetic, and changes through a list, array, or shared reference are
visible to subsequent Python and C++ reads. A C++ consumer converts the current
value to its requested scalar or one-dimensional vector type; incompatible
metadata and out-of-range conversions raise an exception. Python metadata may
use other shapes and containers supported by the HDF5 writer. Objects such as
`None` can be held in memory but have no ALPS HDF5 representation.

Native C++ numeric and Boolean vectors become NumPy arrays when accessed from
Python. Native string vectors become lists so names can be replaced with longer
strings or appended without NumPy's fixed-width string truncation. These
materialized objects retain mutations for subsequent Python and C++ reads.
Explicitly supplied Python lists and NumPy arrays keep their original types.

Integer conversion from text is range checked in the C++ SDK, including when
parameters originate outside Python. Negative text converted to an unsigned
integer now raises an exception instead of wrapping; replace negative textual
sentinels with an explicit value in the target type's range.

The C++ SDK remains independent of Python and nanobind. Python-owned values and
their checkpoint decoder are supplied by the bindings. Rebuild downstream C++
extensions against the SDK from the same source revision as the wheel; the
parameter layout changed during this migration.

Downstream nanobind modules must also use the same nanobind internals ABI as
the installed wheel; otherwise nanobind cannot see pyalps types and aborts the
interpreter at import (`base type "alps::mcbase" not known to nanobind`). The
wheel is built with the nanobind version pinned in `pyproject.toml` and records
its ABI in `pyalps.pyalps_config`; `alps_target_link_pyalps` rejects a
mismatched nanobind at CMake configure time. Install the matching release, for
example `python -m pip install "nanobind==$(python -c 'import
pyalps.pyalps_config as c; print(c.NANOBIND_VERSION)')"`. The compiler's C++
standard library must also match (libc++ on macOS, libstdc++ on Linux).

New HDF5 writes distinguish Boolean and signed-byte values with an
`__alps_type__` attribute while retaining the existing numeric storage format.
Unmarked signed-byte data from old ALPS files retains the legacy Boolean
interpretation. The old format cannot distinguish an unmarked `int8` array
from a Boolean mask; use a typed reader such as h5py when an old dataset is
known to contain signed bytes.

Rectangular mixtures of numeric rows are stored as a single array when every
integer remains exact in the common dtype. If mixing integer widths or mixing
integers with floating-point or complex rows would round a value, the archive
stores the rows separately and reads them back as a list. For example, a
`uint64` row containing `2**63 + 1` alongside an `int64` row retains its exact
integer values instead of silently converting them to `float64`.

`pyalps.mpi` receives Python objects using matched probes, so asynchronous
receives and the wait/test helpers can handle messages larger than mpi4py's
default object receive buffer. This adapter exchanges mpi4py messages;
Boost.MPI's C++ serialization protocol and skeleton/content API are not wire
compatible. Communicating processes must use the same protocol.

Communicator wrappers compare equal when their underlying mpi4py communicators
compare equal. They are intentionally unhashable, matching mpi4py. Unlike the
old Boost.MPI wrappers, they cannot be used as dictionary keys or set members;
applications needing such associations should use explicit application keys.

## Versioning

pyalps does not carry a version of its own. The numeric version is read from
`ALPS_VERSION.txt` at the repository root — the same file
`cmake/ALPSVersion.cmake` reads for `ALPS_VERSION_CORE` — so a release bump is
one edit rather than two that can drift. `test/pyalps/test_wheel_payload.py`
fails if the installed version and that file disagree.

A prerelease label cannot live in that file: `project(VERSION ...)` rejects a
non-numeric version, and neither `find_package()` matching nor the library
SOVERSION has a notion of prerelease ordering. CMake takes it from the
`ALPS_VERSION_PRERELEASE` cache variable; the Python build takes it from the
environment variable of the same name, using the same vocabulary:

| `ALPS_VERSION_PRERELEASE` | version with `ALPS_VERSION.txt` = 2.3.4 |
|---|---|
| unset | `2.3.4` |
| `beta.1` | `2.3.4b1` |
| `alpha.2` | `2.3.4a2` |
| `rc.1` | `2.3.4rc1` |
| `dev.3` | `2.3.4.dev3` |

In GitHub release builds, the provider takes the prerelease label from
`GITHUB_REF` (for example, `refs/tags/v3.0.0-beta.1`). It rejects a tag whose
numeric version differs from `ALPS_VERSION.txt`, or whose label conflicts with
an explicit `ALPS_VERSION_PRERELEASE`. Wheels and source distributions use the
same provider. `python python/pyalps/_build_support/alps_version.py`
prints the version a build would produce, from any working directory.
An sdist preserves its recorded version when rebuilt without the original
build environment.

Note the consequence: because the number is inherited, a Python-only API change
cannot be signalled in the pyalps version alone — it takes a bump of
`ALPS_VERSION.txt`, which moves the whole project.
