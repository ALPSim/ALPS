cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(alpsize NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

# rule for generating single core ising program
add_executable(single ising.cpp single.cpp)
target_link_libraries(single ${ALPS_LIBRARIES})

if (MPI_FOUND)

	# rule for generating mpi ising program
	add_executable(mpi ising.cpp mpi.cpp)
	target_link_libraries(mpi ${ALPS_LIBRARIES})

	# rule for generating mpi pscan ising program
	add_executable(mpi_pscan ising.cpp mpi_pscan.cpp)
	target_link_libraries(mpi_pscan ${ALPS_LIBRARIES})

endif (MPI_FOUND)
