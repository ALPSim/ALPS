#  Copyright Matthias Troyer and Synge Todo 2009 - 2015.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

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
