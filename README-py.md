[![ALPS CI/CD](https://github.com/ALPSim/legacy/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/legacy/actions/workflows/build.yml)

## Python Algorithms and Libraries for Physics Simulations

This is python packages for `Algorithms and Libraries for Physics Simulations` project. For more information check [README.txt](https://pypi.org/project/pyalps/2.3.3/README.txt).

### Installation instruction from binaries

1. pyALPS can be installed on most Linux and MacOS machines from prebuilt binaries available on [PyPI](https://pypi.org/project/pyalps).
pyALPS can be installed using `pip` Python package manager:

```
pip install pyalps
```

Wheels are published for CPython on x86_64 Linux (manylinux) and Apple-silicon macOS; the Python versions covered by a release are listed under "Download files" on PyPI. If no wheel matches your interpreter, `pip` silently falls back to compiling pyALPS **and Boost** from the source distribution, which takes an hour or more and several gigabytes of memory. To make `pip` stop with a clear message instead, install with

```
pip install --only-binary=:all: pyalps
```

### Installation on Windows

pyALPS has no native Windows build. Use the [Windows Subsystem for Linux](https://learn.microsoft.com/windows/wsl/install) (WSL 2) with a Linux distribution such as Ubuntu, create a virtual environment inside WSL with a Python version for which a Linux wheel exists (see above), and run the `pip install` command there.

If the installation is killed while "Building wheel for pyalps", `pip` found no wheel and started a source build that ran out of memory; WSL becoming unresponsive at the same time is the usual symptom. Switch to a supported Python version or, if you do want to build from source inside WSL, first give the VM enough memory by creating `%USERPROFILE%\.wslconfig` on the Windows side:

```
[wsl2]
memory=12GB
swap=32GB
processors=4
```

Run `wsl --shutdown` from PowerShell, reopen the WSL terminal, and follow the source-build instructions below. Note that the pyalps 2.3.3 source distribution on PyPI does not build with scikit-build-core 0.10 or newer (`ERROR: Use build.verbose instead of cmake.verbose`); build from a Git checkout as described below instead.

### Installation instruction from sources

1. Prerequisites
  - CMake >= 3.22
  - Boost sources >= 1.76
  - BLAS/LAPACK
  - HDF5
  - MPI
  - Python >= 3.9
    - Python 3.13 requires Boost version 1.87 or later
    - Earlier versions maybe also work but unsupported
  - C++ compiler with C++17 support (CI covers GCC 11 through 15, Clang 14 through 22, and AppleClang)
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
Boost_SRC_DIR=`pwd`/../boost_1_86_0 python3 -m build --wheel
```
This will download the most recent version of ALPS from the github repository, and build pyALPS python package.

3. Installation

Based on the version of the Python used to build pyALPS, the corresponding Python wheel will be created and stored in `dist` subdirectory. It can be installed using `pip`:
```
pip install dist/pyalps-<specs>.whl
```
