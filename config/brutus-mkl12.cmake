#  Copyright Michele Dolfi 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# We assume the following modules are loaded:
#   module load hdf5 goto2 python/2.7.2 mkl/12.1.2 open_mpi/1.5 cmake
# 
# For the moment python is including goto2 automatically for NumPy, this
# configuration is setting the correct variables to avoid the default CMake FindBLAS.
#
# According to the Intel MKL linker advisor we need:
# * MKL sequential:
#   -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
# * MKL multi-threads:
#   -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm

SET(CMAKE_CXX_COMPILER $ENV{CXX})
SET(CMAKE_C_COMPILER $ENV{CC})
SET(CMAKE_Fortran_COMPILER $ENV{FC})


INCLUDE_DIRECTORIES($ENV{MKLROOT}/include)
SET(LAPACK_LIBRARY "-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -fopenmp -lpthread -lm")
#SET(LAPACK_LIBRARY -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm)

set(HAVE_MKL 1)
SET(BLAS_LIBRARY_INIT 1)
SET(LAPACK_LIBRARY_INIT 1)

