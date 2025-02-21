# - FindMPIForALPS.cmake
# Wrapper to help CMake finding OpenMPI on Mac OS.
# 

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
