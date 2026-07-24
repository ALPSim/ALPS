cmake_minimum_required(VERSION 3.18)
project(ndim_spin NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

set(bench_flags "-march=native -mavx") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${bench_flags}")

# rule for generating the heisenberg example program
add_executable(heisenberg heisenberg.cpp)
target_link_libraries(heisenberg ${ALPS_LIBRARIES})
