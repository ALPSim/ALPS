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

if (ALPS_BUILD_SOURCE)
if(ALPS_PYTHON_WHEEL)

else()
  install(FILES CMakeLists.txt qmc.h DESTINATION alps/applications/qmc)
  install(DIRECTORY checksign looper qwl sse4 worms dwa DESTINATION alps/applications/qmc
          FILES_MATCHING PATTERN "*"  PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE)
endif()
else (ALPS_BUILD_SOURCE)

  add_subdirectory(checksign)
  add_subdirectory(looper)
  add_subdirectory(qwl)
  #add_subdirectory(sse)
  add_subdirectory(sse4)
  add_subdirectory(worms)
  add_subdirectory(dwa)

endif (ALPS_BUILD_SOURCE)
