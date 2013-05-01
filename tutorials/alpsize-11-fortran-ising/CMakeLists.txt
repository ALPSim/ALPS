cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(alpsize NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C, C++, and Fortran compilers
enable_language(C CXX Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")

# rule for generating ising program
add_executable(ising ising.C ising_impl.f90)
target_link_libraries(ising ${ALPS_LIBRARIES} ${ALPS_FORTRAN_LIBRARIES})
