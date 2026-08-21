# pyalps

Python applications and libraries for the Algorithms and Libraries for
Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Install `pyalps[plot]` to use the Matplotlib plotting helpers.
Install `pyalps[mpi]` for the mpi4py-backed `pyalps.mpi` compatibility layer.

The bindings are built as a standalone `scikit-build-core` project using
nanobind. A source build requires Python 3.10 or newer, CMake 3.21 or newer,
Ninja, a C++17 compiler, BLAS/LAPACK, HDF5, and an installed ALPS C++ SDK.
Point `ALPS_DIR` at the SDK's `share/alps` package directory.

The `wheel-deps` CMake preset builds the SDK exactly as the wheel CI does.
From the repository root:

```sh
cmake --preset wheel-deps
cmake --build --preset wheel-deps

ALPS_DIR="$PWD/_build/wheel-deps/install/share/alps" \
  python -m build --wheel bindings/python/pyalps
```

The wheel is written to `bindings/python/pyalps/dist` and can be installed
with `python -m pip install`. With ccache installed, configure with
`cmake --preset wheel-deps -DCMAKE_CXX_COMPILER_LAUNCHER=ccache` and set
`CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"` for the wheel build to
speed up rebuilds.

`PYALPS_BUILD_APPLICATIONS=ON` is the default and preserves the MaxEnt,
CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.

`PYALPS_BUNDLE_APPLICATIONS=ON` is the default and copies the ALPS
application executables (`spinmc`, `dmrg`, `sparsediag`, `loop`, `qwl`, ...)
from the SDK into `pyalps/bin`, together with the SDK's shared libraries in
`pyalps/lib` that their `../lib` RPATH resolves against. `pyalps.tools`
prepends `pyalps/bin` to `PATH`, so this is what makes
`pyalps.runApplication('spinmc', ...)` work from a wheel install — the
`wheel-deps` preset therefore builds the applications. Configure with
`-DPYALPS_BUNDLE_APPLICATIONS=OFF` for a bindings-only wheel; the
`runApplication` helpers then require the executables on `PATH` by other
means.

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

The wheels CI sets it in `[tool.cibuildwheel.environment]`; clear it there to
publish a final release. `python bindings/python/pyalps/_build_support/alps_version.py`
prints the version a build would produce.

Note the consequence: because the number is inherited, a Python-only API change
cannot be signalled in the pyalps version alone — it takes a bump of
`ALPS_VERSION.txt`, which moves the whole project.
