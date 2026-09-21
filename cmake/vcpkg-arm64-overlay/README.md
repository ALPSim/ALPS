# Native Windows ARM64 numerical dependencies

Only the `windows-arm64` preset enables these two packaging adapters. ALPS uses the same imported targets and source code on both Windows architectures.

The pinned upstream vcpkg `lapack` port selects CLAPACK 3.2.1 on Windows ARM64. That C translation uses f2c return conventions, while the accompanying native OpenBLAS uses the modern Fortran ABI. A concrete failure is single-precision Cholesky factorization of `[[4, 2], [2, 3]]`: `SPOTRF` reports success but returns `sqrt(3)` instead of `sqrt(2)` as the second diagonal entry. The ALPS `numerics` test now covers this failure.

The `openblas` adapter installs the official [OpenBLAS 0.3.34 Windows ARM64 LP64 DLL distribution](https://github.com/OpenMathLib/OpenBLAS/releases/tag/v0.3.34), which includes LAPACK 3.12.0. The archive and license are SHA-512 pinned. The upstream CMake imported target is preserved; there is no patched numerical implementation or separate Fortran compiler. The `lapack` adapter exposes that same library through the existing `lapack` CMake target, avoiding CLAPACK.

The upstream package contains a Release C-ABI DLL; it is also the numerical runtime for Debug consumers. Static linkage and non-ARM64 targets are outside this overlay's scope. Ordinary Windows x64 continues to use upstream source ports. Remove these adapters when upstream vcpkg supplies an ABI-compatible modern ARM64 LAPACK provider.
