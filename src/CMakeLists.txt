configure_file(alps/config.h.in ${CMAKE_BINARY_DIR}/src/alps/config.h)
configure_file(ietl/config.h.in ${CMAKE_BINARY_DIR}/src/ietl/config.h)

install(DIRECTORY alps boost ietl mocasito ${CMAKE_BINARY_DIR}/src/alps
        DESTINATION include
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE
	PATTERN ".DS_Store" EXCLUDE
	PATTERN "*.in" EXCLUDE
	PATTERN "Makefile*" EXCLUDE
	)
