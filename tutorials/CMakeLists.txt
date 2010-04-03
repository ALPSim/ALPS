#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

install(DIRECTORY . DESTINATION tutorials COMPONENT tutorials
        FILES_MATCHING PATTERN "*.py" PATTERN "parm*" 
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE 
        PATTERN "*.h5" EXCLUDE PATTERN "*.xml" EXCLUDE 
        PATTERN "*.xsl" EXCLUDE PATTERN "*.run*" EXCLUDE)

