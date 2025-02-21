#  Copyright Matthias Troyer 2009 - 2010.
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

######################################################################
# applications
######################################################################

enable_testing()

if (ALPS_BUILD_SOURCE)
  install(FILES CMakeLists.txt DESTINATION alps/applications)
  install(DIRECTORY mc diag dmft dmrg DESTINATION alps/applications
          FILES_MATCHING PATTERN "*" PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
  add_subdirectory(qmc)
else (ALPS_BUILD_SOURCE)
  add_definitions(${LAPACK_DEFINITIONS})
  add_subdirectory(mc)
  add_subdirectory(qmc)
  add_subdirectory(diag)
  add_subdirectory(dmrg)
  add_subdirectory(dmft)
endif (ALPS_BUILD_SOURCE)


######################################################################
# GUI's
######################################################################

# add_subdirectory(gui/simulator)

