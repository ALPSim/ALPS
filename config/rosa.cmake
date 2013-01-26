#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

SET(CMAKE_CXX_COMPILER CC)
SET(CMAKE_C_COMPILER cc)
SET(CMAKE_Fortran_COMPILER ftn)
SET(NOT_ALPS_BUILD_PYTHON ON)
SET(NOT_BUILD_SHARED_LIBS ON)
SET(BLAS_LIBRARY "/opt/cray/libsci/12.0.00/GNU/47/interlagos/lib/libsci_gnu.a")
SET(LAPACK_LIBRARY "/opt/cray/libsci/12.0.00/GNU/47/interlagos/lib/libsci_gnu.a")
SET(MPI_INCLUDE_PATH "/opt/cray/mpt/default/gni/mpich2-gnu/47/include")
SET(MPI_LIBRARY "/opt/cray/mpt/default/gni/mpich2-gnu/47/lib/libmpich.a")
