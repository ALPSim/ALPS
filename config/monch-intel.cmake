#  Copyright Michele Dolfi 20123.
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
SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -DNDEBUG -DBOOST_DISABLE_ASSERTS" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "" FORCE)

#STRING(REGEX REPLACE "^(.+)/icc$" "\\1" INTEL_BINDIR "${CMAKE_C_COMPILER}")
#MESSAGE(STATUS "Intel binary directory is: ${INTEL_BINDIR}")

SET(CMAKE_AR "xiar" CACHE STRING "" FORCE) 
SET(CMAKE_LINKER "xild" CACHE STRING "" FORCE)

SET(LAPACK_FOUND 1) # with -mkl we have lapack no matter what cmake thinks
SET(HAVE_MKL 1)
SET(BLAS_LIBRARY_INIT 1)
SET(LAPACK_LIBRARY_INIT 1)
SET(BLAS_LIBRARY "-mkl=sequential")
SET(LAPACK_LIBRARY "")
