[![ALPS CI/CD](https://github.com/ALPSim/legacy/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/legacy/actions/workflows/build.yml)

## Python Algorithms and Libraries for Physics Simulations

This is python packages for `Algorithms and Libraries for Physics Simulations` project. For more information check [README.txt](README.txt).

### Installation instruction from binaries

1. pyALPS can be installed on most Linux and MacOS mcachines from prebuilt biniaries available on [PyPi](https://pypi.org/project/pyalps).
pyALPS can be installed using `pip` Python package manager:

```
pip install pyalps
```

### Installation instruction from sources

1. Prerequisites
  - CMake > 3.18
  - Boost sources >= 1.76
  - BLAS/LAPACK
  - HDF5
  - MPI
  - Python >= 3.9
    - Python 3.13 requires Boost version 1.87 or later
    - Earlier versions maybe also work but unsupported
  - C++ compiler (build has been tested on GCC 10.5 through 14.2)
  - GNU Make or Ninja build system

You need to download and unpack boost library:
```
wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
tar -xzf boost_1_86_0.tar.gz
```
Here we download `boost v1.86.0`, we have tested ALPS with versions `1.76.0` and `1.86.0`.

2. Downloading and building sources
```
git clone https://github.com/alpsim/ALPS ALPS
cd ALPS
Boost_ROOT_DIR=`pwd`/../boost_1_86_0 python3 -m build --wheel
```
This will download the most recent version of ALPS from the github repository, and build pyALPS python package.

3. Installation

Based on the version of the Python used to build pyALPS, the corresponding Python wheel will be created and stored in `dist` subdirectory. It can be installed using `pip`:
```
pip install dist/pyalps-<specs>.whl
```
