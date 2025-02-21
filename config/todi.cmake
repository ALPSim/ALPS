#  Copyright Matthias Troyer, Synge Todo and Lukas Gamper 2009 - 2010.
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
