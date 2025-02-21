#  Copyright Michele Dolfi 2012.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
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

