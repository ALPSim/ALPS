#  Copyright Michele Dolfi 2012 - 2013
#            Jan Gukelberger 2012 - 2013
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
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

