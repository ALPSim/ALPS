#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(ising)

enable_language(Fortran)

find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "ALPS version: ${ALPS_VERSION}")
include(${ALPS_USE_FILE})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(ising main.C alps_ising.f90)
target_link_libraries(ising ${ALPS_LIBRARIES} ${ALPS_FORTRAN_LIBRARIES})
