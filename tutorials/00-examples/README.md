# Library examples

These small programs demonstrate individual ALPS APIs. For complete simulation workflows and a recommended learning sequence, start with the [tutorial guide](../README.md).

| Directory | Subject |
| --- | --- |
| [alea](alea/) | Statistical estimates and error analysis |
| [hdf5](hdf5/) | Saving custom C++ types in HDF5 archives |
| [ietl](ietl/) | Iterative eigenvalue methods |
| [model](model/) | Symbolic and numerical model construction |
| [parapack](parapack/) | Simulation workers and parameter scans |
| [sampling](sampling/) | Sampling distributions |
| [scheduler](scheduler/) | Scheduled simulations |
| [fortran](fortran/) | Calling Fortran simulation code through the C++ bridge |

## Build against an installed SDK

These examples require CMake 3.27 or newer; see the [CMake setup instructions](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system provides an older version.

From this directory, build the C++ examples together:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Fortran examples are a separate project and require a Fortran compiler:

```sh
cmake -S fortran -B build-fortran -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps
cmake --build build-fortran --parallel 2
ctest --test-dir build-fortran --output-on-failure
```

Both projects register tests by default; pass `-DALPS_BUILD_TESTING=OFF` to build without registering tests. These standalone builds replace the root `ALPS_BUILD_EXAMPLES` option. To use a local ALPS checkout, build and install its SDK first, then pass that installation prefix in the commands above.

CTest locates the SDK's XML resources. When running an executable directly, use the input files in its example directory; set `ALPS_ROOT` to the SDK prefix if the SDK has been relocated.
