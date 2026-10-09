# Changelog

## Unreleased

- Organize sources, headers, tests and tutorials by module while preserving public include names.

- Export native SDK components and solver libraries through CMake targets. Require CMake 3.27, C++17, external Boost 1.76, HDF5 1.10.5 and LP64 BLAS/LAPACK. MPI defaults to OFF; existing caches retain their configured value. Build tutorials against the installed SDK and use `alps-xml` for XML tools.
