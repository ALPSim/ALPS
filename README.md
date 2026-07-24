[![Build](https://github.com/ALPSim/ALPS/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build.yml)
[![Python wheels](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.txt)

# ALPS — Algorithms and Libraries for Physics Simulations

The ALPS software package aims to provide a set of well tested, robust, and standardized components for numerical simulations of condensed matter systems, including bosonic, fermionic, and spin systems. They consist of a set of components that are used in state-of-the-art high performance codes.

**Project website:** [alps.comp-phys.org](https://alps.comp-phys.org/)

## Installation

### Python

Binary `pyalps` wheels are available for supported Linux and macOS systems:

```sh
python -m pip install pyalps
```

Plotting with `pyalps` requires Matplotlib, which can be installed together with the package:

```sh
python -m pip install "pyalps[plot]"
```

### Build from source

A native build requires CMake 3.18 or newer, a C++14 compiler, HDF5, and BLAS/LAPACK. The bundled CMake presets require CMake 3.21 or newer. MPI is enabled by default when available. The legacy Fortran interface is disabled by default. If a Boost source tree is not supplied, configuration downloads one and therefore requires network access.

Configure a release build with an explicit installation prefix, then build and install it:

```sh
cmake -S . -B _build/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/path/to/alps
cmake --build _build/release --parallel
cmake --install _build/release
```

Alternatively, with CMake 3.21 or newer:

```sh
cmake --preset default
cmake --build --preset default
```

Add `-DALPS_ENABLE_MPI=OFF` to the configure command for a serial-only build. To build the legacy Fortran interface and its examples, add `-DALPS_BUILD_FORTRAN=ON`; this requires a Fortran compiler and the HDF5 Fortran component.

Building the Python bindings from source is a separate step against an installed ALPS C++ SDK; see the [`pyalps` build instructions](bindings/python/pyalps/README.md).

Platform-specific binary, source, and Spack instructions are available on the [ALPS installation website](https://alps.comp-phys.org/documentation/install/).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and contribution guidance.

## License and citations

ALPS is distributed under the terms in [LICENSE.txt](LICENSE.txt). See [CITATION.md](CITATION.md) for citation guidance.
