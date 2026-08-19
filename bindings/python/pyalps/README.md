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
DWA, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.

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
- **Stable ABI (abi3):** not enabled or currently supported. Some binding
  paths still inspect CPython type internals (`tp_name`), and no abi3 build
  runs in CI. Per-version wheels are deliberate; do not add `STABLE_ABI`
  until the code is limited-API clean and CI compiles and imports the
  resulting extensions.
