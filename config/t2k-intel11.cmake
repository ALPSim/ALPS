#  Copyright Haruhiko Matsuo 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See http://www.boost.org/LICENSE_1_0.txt)

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
