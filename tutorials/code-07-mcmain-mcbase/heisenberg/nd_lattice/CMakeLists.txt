cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(heisenberg NONE)

set(ALPS_ROOT_DIR /Users/ricoh/Applications/alps)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

# rule for generating single core ising program
add_executable(single heisenberg.cpp single.cpp)
target_link_libraries(single ${ALPS_LIBRARIES})
