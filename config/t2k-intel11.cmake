#  Copyright Haruhiko Matsuo 2010.
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

# Host:        T2K Tsukuba 
# Compiler:    Intel compiler 11.1
# MPI:         MVAPICH2 v1.4
# BLAS/LAPACK: GENERIC

#
# ALPS Options
#
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(ALPS_BUILD_PYTHON OFF CACHE BOOL ""  FORCE)
set(ALPS_BUILD_APPLICATIONS OFF CACHE BOOL "" FORCE)
set(ALPS_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ALPS_ENABLE_OPENMP ON CACHE BOOL "" FORCE)

#
# Boost 1.42
#
set(Boost_ROOT_DIR $ENV{HOME}/src/boost_1_42_0)

#
# Compiler: Intel 11.1
#
set(CMAKE_CXX_COMPILER /opt/intel/Compiler/11.1/069/bin/intel64/icpc)
set(CMAKE_C_COMPILER /opt/intel/Compiler/11.1/069/bin/intel64/icc)

# compiler flags
add_definitions("-DBOOST_DISABLE_ASSERTS -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX")

# exe linker flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -lrdmacm -libverbs -libumad -lrt" CACHE STRING "" FORCE)

#
# MPI:  MVAPICH2-1
#
set(MPI_COMPILER mpicxx)
set(MPIEXEC_NUMPROC_FLAG "-np")

#
# BLAS and LAPACK: MKL
# 

#
# HDF5
#
set(Hdf5_INCLUDE_DIRS $ENV{HOME}/opt/include)
set(Hdf5_LIBRARY_DIRS $ENV{HOME}/opt/lib)
