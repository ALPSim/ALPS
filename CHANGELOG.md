# Changelog

## Unreleased

- Organize sources, headers, tests and tutorials by module while preserving public include names.

- Export native SDK components and solver libraries through CMake targets. Require CMake 3.27, C++17, external Boost 1.76, HDF5 1.10.5 and LP64 BLAS/LAPACK. MPI defaults to OFF; existing caches retain their configured value. Build tutorials against the installed SDK and use `alps-xml` for XML tools.

- Share the nanobind runtime across Python modules and downstream extensions. Require Python 3.11; use cp312-abi3 for Python 3.12 and newer. Package native SDK libraries, solvers and resources, and expose the `pyalps::runtime` CMake target.

- Select PR validation by affected components and dependency fingerprints. Consolidate installed-artifact checks, MPI and sanitizer jobs; run the broad compatibility matrix on scheduled, manual and release workflows.

- Consolidate native tests with GoogleTest and semantic assertions, retaining explicit serialization fixtures and MPI test registration.
