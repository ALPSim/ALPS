cmake_minimum_required(VERSION 3.0)

project (solver)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
##we also need the blas classes
find_package(Greens REQUIRED)
find_package(Fourier REQUIRED)
find_package(Cluster REQUIRED)

include(${ALPS_USE_FILE})
include_directories(${GREENS_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})

add_library(solver
external_solver.cpp)

target_link_libraries(solver ${GREENS_LIBRARIES} ${FOURIER_LIBRARIES} ${CLUSTER_LIBRARIES})
