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
  if(MPI_FOUND)
	  find_program(MPIEXEC
	    NAMES "openmpirun"
	    # PATHS ${_MPI_PREFIX_PATH}
	    PATH_SUFFIXES bin
	    DOC "Executable for running MPI programs.")
  endif(MPI_FOUND)
  
  if(NOT MPI_FOUND)
    message(STATUS "Forcing MPI compiler to 'mpicxx-openmpi-mp':")
    set(MPI_C_COMPILER mpicc-openmpi-mp)
    set(MPI_CXX_COMPILER mpicxx-openmpi-mp)
    find_package(MPI)
    if(MPI_FOUND)
      find_program(MPIEXEC
        NAMES "mpiexec-openmpi-mp"
        # PATHS ${_MPI_PREFIX_PATH}
        PATH_SUFFIXES bin
        DOC "Executable for running MPI programs.")
    endif(MPI_FOUND)
  endif(NOT MPI_FOUND)
endif(NOT MPI_FOUND AND APPLE)
