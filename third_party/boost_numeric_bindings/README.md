# Boost Numeric Bindings

This is a temporary, trimmed copy of [Numeric Bindings](https://github.com/uBLAS/numeric_bindings), a separately distributed C++ interface to BLAS and LAPACK. It is not part of the installed Boost dependency or ALPS's Python bindings.

The original imported revision was not recorded. Copyright and license notices remain in the individual headers; the accompanying `LICENSE_1_0.txt` contains the Boost Software License.

ALPS uses these headers in its public numerical templates, IETL, and solvers. `ALPS::headers` supplies the build include directory, and the SDK installs the remaining headers under `include/boost/numeric/bindings`. Removing them requires migrating those consumers first.

Local changes:

- The Fortran integer types have explicit 32-bit and 64-bit widths. ALPS requires the 32-bit LP64 interface with lowercase underscore symbols.
- Unreferenced headers have been removed, including the Eigen, GLAS, MTL, UMFPACK, Boost.Array, and Boost.MultiArray adapters. The retained include graph covers conditional backend includes as well as direct includes.

The intended replacement is a small ALPS numerical interface backed by supported BLAS/LAPACK APIs. Keep this directory isolated until that migration is complete; avoid adding new consumers of the vendored API.
