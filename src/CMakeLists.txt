#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

configure_file(alps/config.h.in ${CMAKE_BINARY_DIR}/src/alps/config.h)
configure_file(ietl/config.h.in ${CMAKE_BINARY_DIR}/src/ietl/config.h)

install(DIRECTORY alps boost ietl mocasito ${CMAKE_BINARY_DIR}/src/alps
        DESTINATION include
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE
	PATTERN ".DS_Store" EXCLUDE
	PATTERN "*.in" EXCLUDE
	PATTERN "Makefile*" EXCLUDE
    COMPONENT headers
	)
