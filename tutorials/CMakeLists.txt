#  Copyright Matthias Troyer and Synge Todo 2009 - 2015.
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

if(NOT ALPS_PYTHON_WHEEL AND NOT ALPS_BUILD_LIBS_ONLY)
if (ALPS_BUILD_SOURCE)
  if (ALPS_BUILD_TUTORIALS_PACKAGE)
    set (destdir .)
  else (ALPS_BUILD_TUTORIALS_PACKAGE)
    set(destdir alps/tutorials)
  endif (ALPS_BUILD_TUTORIALS_PACKAGE)
    
  install(FILES CMakeLists.txt cleanup.sh DESTINATION ${destdir})

  install(DIRECTORY .
          DESTINATION ${destdir}
          FILES_MATCHING PATTERN "*.py" PATTERN "*.ipynb" PATTERN "*.sh" PATTERN "parm*" PATTERN "*params" PATTERN "*.ip" PATTERN "*.op" PATTERN "*.dat"
          PATTERN "*.parm" PATTERN "*.vt" PATTERN "*.xml" PATTERN "*input*" PATTERN "*.pvsm"
	  PATTERN "*.h" PATTERN "*.c" PATTERN "*.cpp" PATTERN "*.C" PATTERN "*.f" PATTERN "*.f90" PATTERN "Makefile" PATTERN "CMakeLists.txt" 
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE PATTERN "*.chkp" EXCLUDE
        PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE PATTERN "*clone*" EXCLUDE PATTERN "dump.xml" EXCLUDE)
else (ALPS_BUILD_SOURCE)
  install(DIRECTORY . DESTINATION tutorials COMPONENT tutorials
          FILES_MATCHING PATTERN "*.py" PATTERN "*.ipynb" PATTERN "*.sh" PATTERN "parm*" PATTERN "*params" PATTERN "*.ip" PATTERN "*.op" PATTERN "*.dat"
          PATTERN "*.parm" PATTERN "*.vt" PATTERN "*.xml" PATTERN "*input*" PATTERN "*.pvsm"
	  PATTERN "*.h" PATTERN "*.c" PATTERN "*.cpp" PATTERN "*.C" PATTERN "*.f" PATTERN "*.f90" PATTERN "Makefile" PATTERN "CMakeLists.txt"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE PATTERN "*.chkp" EXCLUDE
          PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE PATTERN "*clone*" EXCLUDE PATTERN "dump.xml" EXCLUDE)
endif (ALPS_BUILD_SOURCE)
endif()