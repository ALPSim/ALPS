# pyalps

Legacy-compatible Python bindings for ALPS, built as a standalone
`scikit-build-core` project using nanobind. The C++ ALPS library must be
built and installed separately; point `ALPS_DIR` at its `share/alps`
package directory when building this wheel.

From the repository root:

```sh
cmake -S . -B _build/alps -G Ninja \
  -DCMAKE_INSTALL_PREFIX="$PWD/_build/install" \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DALPS_BUILD_LIBS_ONLY=ON
cmake --build _build/alps --target install

ALPS_DIR="$PWD/_build/install/share/alps" \
  CMAKE_ARGS="-DCMAKE_CXX_COMPILER_LAUNCHER=ccache" \
  python -m build --wheel
```

`PYALPS_BUILD_APPLICATIONS=ON` is the default and preserves the MaxEnt,
DWA, CT-HYB, and CT-INT extension modules. Set it to `OFF` through CMake
configuration for a smaller core-only developer build.
