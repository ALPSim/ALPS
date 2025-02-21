#  Copyright Jan Gukelberger 2012.
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
#
# We assume the following modules are loaded:
#   module load hdf5 goto2 python/2.7.2 intel/12 open_mpi/1.[56] cmake
# Do NOT load module mkl, it is already included with intel/12.
#
# icpc option "-fast" enables "-static" which does not work well with cmake's
# library search when there are dynamic libraries around. 
# Therefore, in CMAKE_CXX_FLAGS_RELEASE we manually set all other options
# included in "-fast". 

SET(CMAKE_CXX_COMPILER $ENV{CXX})
SET(CMAKE_C_COMPILER $ENV{CC})
SET(CMAKE_Fortran_COMPILER $ENV{FC})
SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -mkl -DNDEBUG -DBOOST_DISABLE_ASSERTS" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -mkl" CACHE STRING "" FORCE)

STRING(REGEX REPLACE "^(.+)/icc$" "\\1" INTEL_BINDIR "${CMAKE_C_COMPILER}")
#MESSAGE(STATUS "Intel binary directory is: ${INTEL_BINDIR}")

SET(CMAKE_AR "${INTEL_BINDIR}/xiar" CACHE STRING "" FORCE) 
SET(CMAKE_LINKER "${INTEL_BINDIR}/xild" CACHE STRING "" FORCE)

SET(LAPACK_FOUND 1) # with -mkl we have lapack no matter what cmake thinks
SET(HAVE_MKL 1)
SET(BLAS_LIBRARY_INIT 1)
SET(LAPACK_LIBRARY_INIT 1)
SET(BLAS_LIBRARY "")
SET(LAPACK_LIBRARY "")
