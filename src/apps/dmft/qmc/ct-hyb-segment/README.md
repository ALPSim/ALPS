#### CT-HYB-SEGMENT
`CT-HYB-SEGMENT` is the hybridization expansion continuous-time quantum Monte Carlo code, implemented and described in:
H. Hafermann, P. Werner, E. Gull, CPC 184, 1280 (2013).
This version of the code is its adaptation to *ALPSCore* library. 

#### Dependencies
1. ALPSCore (http://alpscore.org)
2. LAPACK
3. NFFT3 (https://www-user.tu-chemnitz.de/~potts/nfft/)

#### Usage 
```
alps_cthyb parameters
```
See `alps_cthyb --help` for the list of parameters

Python bindings are temporarily disabled

#### Installation

Replace the variables with `${}` to the proper values for your system.

```ShellSession
$ git clone https://github.com/ALPSCore/CT-HYB-SEGMENT.git
$ cd CT-HYB-SEGMENT
$ CTHYB_INSTALL_DIR=${where-to-install}
$ mkdir build
$ cd build
$ NFFT3_DIR=${NFFT3-install-dir} CC=${C-compiler} CXX=${C++-compiler} cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DALPSCore_DIR=${ALPSCore_INSTALL_DIR}/share/ALPSCore \
    -DCMAKE_INSTALL_PREFIX=$CTHYB_INSTALL_DIR
$ make
$ make install
```
