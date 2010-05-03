#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_SOURCE)
  install(FILES CMakeLists.txt cleanup.sh DESTINATION alps/tutorials)

  install(DIRECTORY 
            00-introduction
            01-autocorrelations
            02-mc-qmc-susceptibilities
            03-qmc-magnetization
            04-qmc-measurements
            05-qmc-bosons
            06-qwl
            07-sparsediag
            08-fulldiag
            09-particle-in-box-dmrg
            10-dmrg
            mc-01-autocorrelations
            mc-02-susceptibilities
          DESTINATION alps/tutorials
          FILES_MATCHING PATTERN "*.py" PATTERN "parm*" PATTERN "*.vt" PATTERN "*.xml"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE
          PATTERN "*clone*" EXCLUDE
          PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE)
else (ALPS_BUILD_SOURCE)
  install(DIRECTORY . DESTINATION tutorials COMPONENT tutorials
          FILES_MATCHING PATTERN "*.py" PATTERN "parm*" PATTERN "*.vt" PATTERN "*.xml"
          PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
          PATTERN "*.h5" EXCLUDE  EXCLUDE  PATTERN "parm*xml" EXCLUDE
          PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE)
endif (ALPS_BUILD_SOURCE)

