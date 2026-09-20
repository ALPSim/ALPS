The `vcpkg-gfortran` overlay is based on vcpkg baseline
`e6f9e70a29a3e80a1fc510d8503304315447112f` (port version `3#3`).

Only the `windows-x64-on-arm64` preset enables it. The ordinary `windows-x64`
preset and GitHub Windows jobs use upstream ports without an overlay.

It removes the host-architecture check and its unused `MINGW_W`/`MSYS_HOST`
variables. The port copies runtime DLLs using `VCPKG_TARGET_ARCHITECTURE`;
the host architecture does not affect those paths. Windows ARM64 can execute
the downloaded x64 MinGW compiler under Windows emulation, allowing the x64
LAPACK dependency to build on this development host.

This does not provide an ARM64 Fortran compiler or native ARM64 ALPS binaries.
Windows 11 ARM64 runs the x64 compiler and applications under emulation.
For native ARM64 output, use the `windows-arm64` preset and its separate
[numerical dependency adapter](../vcpkg-arm64-overlay/README.md).
Remove the overlay when the pinned upstream port drops this check.

Upstream context: [vcpkg issue #45721](https://github.com/microsoft/vcpkg/issues/45721).
Overlay ports are an [official vcpkg extension mechanism](https://learn.microsoft.com/en-us/vcpkg/concepts/overlay-ports).
