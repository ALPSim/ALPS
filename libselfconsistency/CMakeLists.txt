cmake_minimum_required(VERSION 3.0)

project (selfconsistency)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
#we also need the blas classes
find_package(BlasClasses REQUIRED)
find_package(Greens REQUIRED)
find_package(Fourier REQUIRED)
find_package(Cluster REQUIRED)
find_package(Solver REQUIRED)

include(${ALPS_USE_FILE})
include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GREENS_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${SOLVER_INCLUDE_DIR})

add_library (selfconsistency
hilbert_transform.cpp        selfconsistency_CDMFT.cpp    selfconsistency_DCAPlus.cpp
selfconsistency.cpp          selfconsistency_DCA.cpp
)

target_link_libraries(selfconsistency ${SOLVER_LIBRARIES} ${FOURIER_LIBRARIES} ${GREENS_LIBRARIES} ${CLUSTER_LIBRARIES})
