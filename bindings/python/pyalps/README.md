# pyalps

Python applications and libraries for the Algorithms and Libraries for
Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Install `pyalps[plot]` to use the Matplotlib plotting helpers.

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
