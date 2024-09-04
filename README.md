[![ALPS CI/CD](https://github.com/ALPSim/legacy/actions/workflows/build.yml/badge.svg)](https://github.com/ALPSim/legacy/actions/workflows/build.yml)

## Algorithms and Libraries for Physics Simulations

This is the legacy reposiory. For more information check [README.txt](README.txt).

### Installation instruction

1. Prerequisites
  - CMake > 2.8
  - Boost sources >= 1.68
  - BLAS/LAPACK
  - HDF5
  - MPI
  - Python >= 3.9 (earlier versions maybe also work but unsupported)
  - C++ compiler (build has been tested on GCC 10.5, GCC 11.4, and GCC 12.3)
  - GNU Make or Ninja build system

You need to download and unpack boost library:
```
wget https://archives.boost.io/release/1.76.0/source/boost_1_76_0.tar.gz
tar -xzf boost_1_76_0.tar.gz
```
Here we download `boost v1.76.0`, we have tested ALPS with versions `1.76.0` and `1.81.0`.

2. Downloading and building sources
```
git clone https://github.com/alpsim/legacy alps_legacy
cmake -S alps_legacy -B alps_build -DCMAKE_INSTALL_PREFIX=</path/to/install/dir> \
      -DBoost_ROOT_DIR=`pwd`/boost_1_76_0                                        \
      -DCMAKE_CXX_FLAGS="-std=c++14 -fpermissive"
cmake --build alps_build -j 8
cmake --build alps_build -t test
```
This will download the most recent version of ALPS from the github repository, build it, and run unit tests.

3. Installation
```
cmake --build alps_build -t install
```
