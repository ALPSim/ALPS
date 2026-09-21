[![Build](https://github.com/ALPSim/ALPS/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build.yml) [![Python wheels](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/build_wheels.yml) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.txt)

# ALPS — Algorithms and Libraries for Physics Simulations

The ALPS software package aims to provide a set of well tested, robust, and standardized components for numerical simulations of condensed matter systems, including bosonic, fermionic, and spin systems. They consist of a set of components that are used in state-of-the-art high performance codes.

**Project website:** [alps.comp-phys.org](https://alps.comp-phys.org/)

## Learn ALPS

Start with the [tutorial guide](tutorials/README.md): run your first simulation, choose a numerical method, or learn to develop with the ALPS libraries. It includes recommended learning paths, notebooks, and focused API examples.

## Installation

For current binary, source, and Spack installation instructions, see the [ALPS installation website](https://alps.comp-phys.org/install/).

Source builds require **CMake 4.3 or newer**. If your system provides an older version, follow the [CMake and Ninja setup](CONTRIBUTING.md#install-cmake-and-ninja) for pip installation on Linux, macOS, and Windows, or official binary downloads. For the full build instructions, including MSVC/vcpkg presets for Windows x64 and ARM64 and the exported CMake SDK target, see [the development setup](CONTRIBUTING.md#getting-started-with-the-code).

Python projects live under [`python/`](python/README.md); see the [pyalps build instructions](python/pyalps/README.md) to build the current package against an installed C++ SDK.

See [CHANGELOG.md](CHANGELOG.md) for release notes and upgrade guidance.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and contribution guidance.

## License and citations

ALPS is distributed under the terms in [LICENSE.txt](LICENSE.txt). See [CITATION.md](CITATION.md) for citation guidance.
