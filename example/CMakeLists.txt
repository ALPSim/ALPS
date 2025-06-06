#  Copyright Matthias Troyer and Synge Todo 2010-2011.
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


if (ALPS_BUILD_SOURCE)

  INSTALL(FILES CMakeLists.txt DESTINATION alps/example)  
  install(DIRECTORY alea parapack scheduler sampling model ietl fortran DESTINATION alps/example
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)

else (ALPS_BUILD_SOURCE)

  enable_testing()

  add_subdirectory(alea)
  add_subdirectory(parapack)
  add_subdirectory(scheduler)
  add_subdirectory(ngs)
  add_subdirectory(sampling)
  add_subdirectory(model)
  add_subdirectory(ietl)
  add_subdirectory(hdf5)
  if (ALPS_BUILD_FORTRAN)
    add_subdirectory(fortran)
  endif (ALPS_BUILD_FORTRAN)

endif (ALPS_BUILD_SOURCE)
