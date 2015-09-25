cmake_minimum_required(VERSION 2.8)

project (DMFT)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package(LAPACK REQUIRED)
find_package(BlasClasses REQUIRED)
find_package(GFTools REQUIRED)
#find_package(Fourier REQUIRED)
find_package(Cluster REQUIRED)
find_package(Solver REQUIRED)
find_package(Selfconsistency REQUIRED)

include(UseGtest)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
#include_directories(${GREENS_INCLUDE_DIR})
#include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${SOLVER_INCLUDE_DIR})
include_directories(${SELFCONSISTENCY_INCLUDE_DIR})

add_executable(dmft
main.cpp
)

target_link_libraries(dmft ${SELFCONSISTENCY_LIBRARIES} ${SOLVER_LIBRARIES} ${CLUSTER_LIBRARIES} #${FOURIER_LIBRARIES} 
#${GREENS_LIBRARIES} 
${CLUSTER_LIBRARIES} ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS dmft RUNTIME DESTINATION "bin/")
install(FILES cluster1d.xml cluster3d.xml cluster.xml DESTINATION "share/")
