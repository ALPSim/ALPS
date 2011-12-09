#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

SET(CMAKE_CXX_COMPILER CC)
SET(CMAKE_C_COMPILER cc)
SET(CMAKE_Fortran_COMPILER ftn)
SET(NOT_ALPS_BUILD_PYTHON ON)
SET(NOT_BUILD_SHARED_LIBS ON)
SET(BLAS_LIBRARY "/opt/xt-libsci/default/cray/74/interlagos/lib/libsci_cray.a")
SET(LAPACK_LIBRARY "/opt/xt-libsci/default/cray/74/interlagos/lib/libsci_cray.a")
SET(MPI_INCLUDE_PATH "/opt/cray/mpt/5.3.4/xt/gemini/mpich2-cray/73/include")
SET(MPI_Fortran_INCLUDE_PATH "/opt/cray/mpt/5.3.4/xt/gemini/mpich2-cray/73/include")
SET(MPI_LIBRARY "/opt/cray/mpt/5.3.4/xt/gemini/mpich2-cray/73/lib/libmpich_cray.a")
SET(MPI_Fortran_LIBRARIES "/opt/cray/mpt/5.3.4/xt/gemini/mpich2-cray/73/lib/libmpich_cray.a")
SET(HDF5_INCLUDE_DIR "/opt/cray/hdf5/default/hdf5-cce/include/")
SET(HDF5_LIBRARY "/opt/cray/hdf5/default/hdf5-cce/lib/libhdf5.a")
