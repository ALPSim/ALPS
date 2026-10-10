[![Build](https://github.com/ALPSim/ALPS/actions/workflows/ci.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/ci.yml) [![Compatibility](https://github.com/ALPSim/ALPS/actions/workflows/compatibility.yml/badge.svg)](https://github.com/ALPSim/ALPS/actions/workflows/compatibility.yml) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.txt)

# ALPS — Algorithms and Libraries for Physics Simulations

The ALPS software package aims to provide a set of well tested, robust, and standardized components for numerical simulations of condensed matter systems, including bosonic, fermionic, and spin systems. They consist of a set of components that are used in state-of-the-art high performance codes.

**Project website:** [alps.comp-phys.org](https://alps.comp-phys.org/)

## Learn ALPS

Start with the [tutorial guide](tutorials/README.md): run your first simulation, choose a numerical method, or learn to develop with the ALPS libraries. It includes recommended learning paths, notebooks, and focused API examples.

## Installation

Build this checkout with standard CMake and Python packaging commands; see [CONTRIBUTING.md](CONTRIBUTING.md#getting-started-with-the-code) for dependencies and platform instructions. Use your existing dependency provider and toolchain. The [ALPS installation website](https://alps.comp-phys.org/install/) also covers upstream releases and Spack; its combined Boost.Python source-build instructions use a different interface from this SDK checkout.

Source builds require **CMake 3.27 or newer**. If your system provides an older version, follow the [CMake and Ninja setup](CONTRIBUTING.md#install-cmake-and-ninja) for pip installation or official binary downloads. For the full build instructions and exported CMake SDK target, see [the development setup](CONTRIBUTING.md#getting-started-with-the-code).

The [pyalps build instructions](python/pyalps/README.md) explain how to build the Python package against an installed C++ SDK.

See [CHANGELOG.md](CHANGELOG.md) for release notes and upgrade guidance.

## Develop from source

After installing the [prerequisites](CONTRIBUTING.md#prerequisites), configure, build and test the SDK:

```sh
cmake --preset default
cmake --build --preset default --parallel 2
ctest --preset default
cmake --install _build/default
```

Then follow the [editable Python installation](python/pyalps/README.md#editable-development), pointing `ALPS_DIR` at `_build/default/install/share/alps`.

The [module layout](src/alps/README.md) prepares MaxEnt, HDF5 and typed params for ALPSCore reconciliation. Utilities, HDF5, params, Osiris, XML and command-line parsing have separate exported libraries; both the MaxEnt solver and executable link components without `ALPS::alps`. Public include names and scientific algorithms are preserved. `ALPS::headers` still provides a shared compile interface, so source ownership does not imply that every module is independent.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development and contribution guidance.

## License and citations

ALPS is distributed under the terms in [LICENSE.txt](LICENSE.txt). See [CITATION.md](CITATION.md) for citation guidance, generated from the bibliography in [CITATION.cff](CITATION.cff) and application policy in [CITATIONS.yaml](CITATIONS.yaml). Simulation applications print the relevant guidance at startup; use `--citations` to display it without input files. `--help` and `--license` show usage and legal terms separately.
