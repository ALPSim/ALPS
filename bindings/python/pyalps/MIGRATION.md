# Migrating from the Boost.Python pyalps build

The public Python API is preserved wherever it maps to native ALPS values.
The nanobind build intentionally does not keep arbitrary Python objects inside
`alps::params` or the C++ library.

## Parameters

`pyalps.ngs.params` accepts native booleans, 32-bit integers, floating-point
and complex numbers, strings, homogeneous Python sequences, NumPy scalar
arrays, and one-dimensional NumPy arrays. Values are copied into native C++
storage. Multidimensional arrays, `None`, dictionaries, arbitrary objects, and
integers outside ALPS' 32-bit parameter range raise `TypeError` rather than
being retained as opaque Python objects. Sequence values are returned as
lists, irrespective of whether the input was a list, tuple, or NumPy array.

## MPI

Install `pyalps[mpi]` to use `pyalps.mpi`. The module provides the commonly
used Boost.MPI Python surface (`world`, `rank`, `size`, `Communicator`,
point-to-point methods, collectives, status/request names, and `Timer`) on top
of mpi4py. Ordinary pyalps wheels remain independent of any MPI runtime.

The historical `mcbase(..., communicator)` argument is still accepted. It is
ignored, as it was by the Boost.Python wrapper; `alps::mcbase` itself has no
communicator constructor. Use `pyalps.mpi` for Python communication and ALPS'
C++ MPI adapters for MPI-aware C++ simulations.

Boost.MPI's Python-object serialization bridge and skeleton/content API are
not reproduced. Hybrid applications should use mpi4py's typed buffer API or
an application-specific native C++ protocol.

## Compiled module paths and DWA vectors

Legacy paths such as `pyalps.pyalea_c` and `pyalps.dwa_c` remain aliases of
the extensions now stored under `pyalps._ext`. The preferred stable import is
still `pyalps.cxx.pyalea_c` for core extensions and `pyalps.dwa` for DWA.

DWA's former `std_vector_*` constructors are compatibility aliases for
Python's `list`. DWA methods return list snapshots, which avoids exposing
mutable C++ container proxies and accepts ordinary Python sequences directly.

## Exporting downstream C++ simulations

The public header `<alps/ngs/detail/export_sim_to_python.hpp>` now implements
the export helper with nanobind while keeping the
`ALPS_EXPORT_SIM_TO_PYTHON` macro. Change the module declaration in an old
export source from:

```cpp
BOOST_PYTHON_MODULE(my_sim) {
```

to:

```cpp
#include <nanobind/nanobind.h>
NB_MODULE(my_sim, m) {
```

and keep the existing export macro call. See
`tutorials/ngs/5_export_python` for a complete standalone CMake build.

After creating the nanobind target, link it with the installed SDK helper:

```cmake
include("${ALPS_PYTHON_USE_FILE}")
alps_target_link_pyalps(my_sim PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
```

Do not link a wheel consumer directly to a second system `libalps`/HDF5
stack. Repaired wheels carry private shared libraries, and stateful values
such as HDF5 handles are valid only in the library image that created them.
The helper selects the wheel's exact runtime when present, retains normal SDK
linking for source installs, and makes direct `import my_sim` work on macOS.

The removed `alps/python/numpy_array.hpp` API should be replaced with
`nanobind::ndarray` or nanobind's STL casters. The old
`alps/hdf5/python.hpp` operators accepted `boost::python::object` and have no
Python-object-free equivalent; use typed `alps::hdf5::archive` operations in
C++ or `pyalps.hdf5` at the Python boundary.
