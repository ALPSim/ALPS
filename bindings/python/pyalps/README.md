# pyalps

Python applications and libraries for the Algorithms and Libraries for
Physics Simulations (ALPS) project. Binary wheels are available from PyPI:

```sh
python -m pip install pyalps
```

Install `pyalps[plot]` to use the Matplotlib plotting helpers.

The bindings are built as a standalone `scikit-build-core` project using
nanobind. A source build requires Python 3.10 or newer, CMake 3.18 or newer,
a C++17 compiler, BLAS/LAPACK, HDF5, and an installed ALPS C++ SDK. Point
`ALPS_DIR` at the SDK's `share/alps` package directory.

From the repository root:

```sh
cmake -S . -B _build/alps -G Ninja \
  -DCMAKE_INSTALL_PREFIX="$PWD/_build/install" \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DALPS_ENABLE_MPI=OFF \
  -DALPS_BUILD_LIBS_ONLY=ON
cmake --build _build/alps --target install

ALPS_DIR="$PWD/_build/install/share/alps" \
  CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache" \
  python -m build --wheel bindings/python/pyalps
```

The wheel is written to `bindings/python/pyalps/dist` and can be installed
with `python -m pip install`.

`PYALPS_BUILD_APPLICATIONS=ON` is the default and preserves the MaxEnt,
DWA, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.
