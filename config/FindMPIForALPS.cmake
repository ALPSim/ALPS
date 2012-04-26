# - FindMPIForALPS.cmake
# Wrapper to help CMake finding OpenMPI on Mac OS.
# 

#  Copyright Michele Dolfi 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


# First try with standard configuration
find_package(MPI)

if(NOT MPI_FOUND AND APPLE)
  message(STATUS "Forcing MPI compiler to 'openmpicxx':")
  set(MPI_C_COMPILER openmpicc)
  set(MPI_CXX_COMPILER openmpicxx)
  find_package(MPI)
endif(NOT MPI_FOUND AND APPLE)
