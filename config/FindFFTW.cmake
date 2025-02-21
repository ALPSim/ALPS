#  Copyright Olivier Parcollet 2010.
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
# This module looks for fftw.
# It sets up : FFTW_INCLUDE_DIR, FFTW_LIBRARIES
# 

SET(TRIAL_PATHS
 /usr/include
 /usr/local/include
 /opt/local/include
 /sw/include
 )
FIND_PATH(FFTW_INCLUDE_DIR fftw3.h ${TRIAL_PATHS} DOC "Include for FFTW")

SET(TRIAL_LIBRARY_PATHS
 /usr/lib 
 /usr/local/lib
 /opt/local/lib
 /sw/lib
 )

SET(FFTW_LIBRARIES "FFTW_LIBRARIES-NOTFOUND" CACHE STRING "FFTW library")
# Try to detect the lib
FIND_LIBRARY(FFTW_LIBRARIES fftw3 ${TRIAL_LIBRARY_PATHS} DOC "FFTW library")

mark_as_advanced(FFTW_INCLUDE_DIR)
mark_as_advanced(FFTW_LIBRARIES)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDE_DIR)

