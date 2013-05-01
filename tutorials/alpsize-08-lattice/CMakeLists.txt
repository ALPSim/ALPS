cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(alpsize NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

# rule for generating lattice program
add_executable(lattice lattice.C)
target_link_libraries(lattice ${ALPS_LIBRARIES})
add_alps_test(lattice)

# rule for generating wolff program
add_executable(wolff wolff.C)
target_link_libraries(wolff ${ALPS_LIBRARIES})
add_alps_test(wolff)
