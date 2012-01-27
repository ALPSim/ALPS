#  Copyright Matthias Troyer 2009 - 2010.
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
          FILES_MATCHING PATTERN "*.py" PATTERN "parm*" PATTERN "*.vt" PATTERN "*.xml" PATTERN "*input*" PATTERN "*parm"
	  PATTERN "*.cpp" PATTERN "Makefile"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE
          PATTERN "*clone*" EXCLUDE
          PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE)
else (ALPS_BUILD_SOURCE)
  install(DIRECTORY . DESTINATION tutorials COMPONENT tutorials
          FILES_MATCHING PATTERN "*.py" PATTERN "parm*" PATTERN
          "*.parm" PATTERN "*.vt" PATTERN "*.xml" PATTERN "*input*"
	  PATTERN "*.cpp" PATTERN "Makefile"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE
          PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE PATTERN
          "dump.xml" EXCLUDE)
endif (ALPS_BUILD_SOURCE)

