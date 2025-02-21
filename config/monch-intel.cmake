#  Copyright Michele Dolfi 2012 - 2013
#            Jan Gukelberger 2012 - 2013
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
#   module load gcc/4.7.3 intel python hdf5 cmake
# we load gcc 4.7.3 to have updated stdlib
#
# icpc option "-fast" enables "-static" which does not work well with cmake's
# library search when there are dynamic libraries around. 
# Therefore, in CMAKE_CXX_FLAGS_RELEASE we manually set all other options
# included in "-fast". 

SET(CMAKE_CXX_COMPILER icpc)
SET(CMAKE_C_COMPILER icc)
SET(CMAKE_Fortran_COMPILER ifort)
SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -DNDEBUG -DBOOST_DISABLE_ASSERTS" CACHE STRING "Flags used by the compiler during release builds")
#STRING(REGEX REPLACE "^(.+)/icc$" "\\1" INTEL_BINDIR "${CMAKE_C_COMPILER}")
#MESSAGE(STATUS "Intel binary directory is: ${INTEL_BINDIR}")

SET(CMAKE_AR "xiar" CACHE PATH "Path to a program.") 
SET(CMAKE_LINKER "xild" CACHE PATH "Path to a program.")

SET(LAPACK_FOUND 1) # with -mkl we have lapack no matter what cmake thinks
SET(HAVE_MKL 1)
SET(BLAS_LIBRARY_INIT 1)
SET(LAPACK_LIBRARY_INIT 1)
IF(ALPS_USE_MKL_PARALLEL)
    SET(CMAKE_CXX_FLAGS "-mkl=parallel" CACHE STRING "Flags used by the compiler during all build types")
    SET(BLAS_LIBRARY "-mkl=parallel")
    SET(ALPS_ENABLE_OPENMP ON CACHE BOOL "Enable OpenMP parallelization") # parallel mkl needs -openmp
ELSE(ALPS_USE_MKL_PARALLEL)
    SET(CMAKE_CXX_FLAGS "-mkl=sequential" CACHE STRING "Flags used by the compiler during all build types")
    SET(BLAS_LIBRARY "-mkl=sequential")
ENDIF(ALPS_USE_MKL_PARALLEL)
SET(LAPACK_LIBRARY "")

