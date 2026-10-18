# Boost Numeric Bindings

This is ALPS's vendored subset of [Numeric Bindings](https://github.com/uBLAS/numeric_bindings), a separately distributed C++ interface to BLAS and LAPACK. It is not part of the installed Boost dependency or ALPS's Python bindings. The headers are included in the repository so building ALPS does not require fetching or installing another dependency.

## Provenance and local changes

The upstream revision of the original import is unknown. ALPS commit [`02149a3a5416951d5439ef4dc7b57dc5290fba41`](https://github.com/ALPSim/ALPS/commit/02149a3a5416951d5439ef4dc7b57dc5290fba41) (November 1, 2010, “Added current bindings snapshot”) placed the snapshot in `src/boost/numeric/bindings`. ALPS subsequently modified it. The baseline for this build modernization is the `bindings/boost/numeric/bindings` tree at ALPS commit [`c22bfd7010511826350523418fb11c6806a87d6a`](https://github.com/ALPSim/ALPS/commit/c22bfd7010511826350523418fb11c6806a87d6a); it should not be identified as an unmodified upstream release.

Changes from that ALPS baseline:

- Relocate the headers to this directory and supply a copy of the Boost Software License.
- Give the Fortran integer and logical types explicit 32-bit and 64-bit widths in `detail/config/fortran.hpp`. ALPS requires the 32-bit LP64 interface with lowercase underscore symbols.
- Remove 295 headers outside ALPS's transitive include graph, including unused adapters and LAPACK wrappers. ALPS includes specific LAPACK routines or the existing Fortran declarations instead of the LAPACK umbrella headers. Retain includes in every conditional branch, including alternative backend declarations; the BLAS umbrella remains unchanged.

Copyright and license notices remain in the individual headers; the accompanying `LICENSE_1_0.txt` contains the Boost Software License and is installed with the SDK headers.

## Maintenance

ALPS uses these headers in its public numerical templates, IETL, and solvers. `ALPS::headers` supplies the build include directory, and the SDK installs the remaining headers under `include/boost/numeric/bindings`. Removing them requires migrating those consumers first.

After changing this subset, run `python -m pytest tests/ci/test_numeric_bindings.py` from the repository root, rebuild the SDK, and run its numerical tests and installed-SDK consumer checks. The include audit checks literal includes across conditional branches; it does not establish that every alternative backend is supported. Document any further patches or snapshot replacement here. Keep the snapshot in the repository rather than downloading an unpinned version during configuration.

A possible future replacement is an ALPS numerical interface backed by supported BLAS/LAPACK APIs. Its backend and dependency choices need a separate evaluation and numerical validation. Until then, keep this directory isolated and avoid adding new consumers of the vendored API.
