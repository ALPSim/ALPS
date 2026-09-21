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

These examples require CMake 4.3 or newer; see the [CMake setup instructions](../../CONTRIBUTING.md#install-cmake-and-ninja) if your system provides an older version.

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

CTest locates the SDK's XML resources and Windows DLLs. When running an executable directly, use the input files in its example directory; set `ALPS_ROOT` to the SDK prefix if the SDK has been relocated. On Windows, also add the SDK's `bin` directory to `PATH`.

## Build with the ALPS source tree

Configure the repository with `-DALPS_BUILD_EXAMPLES=ON` to include the C++ examples in its build and test suite. Fortran examples use the standalone build above.
