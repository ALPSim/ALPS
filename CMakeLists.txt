cmake_minimum_required(VERSION 2.8)

project (DMFT)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
find_package(LAPACK REQUIRED)
find_package(BlasClasses REQUIRED)
find_package(Greens REQUIRED)
find_package(Fourier REQUIRED)
find_package(Cluster REQUIRED)

include(${ALPS_USE_FILE})
include(UseGtest)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GREENS_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})

add_executable(dmft
hilbert_transform.cpp        selfconsistency_CDMFT.cpp    selfconsistency_DCAPlus.cpp
external_solver.cpp          selfconsistency.cpp          selfconsistency_DCA.cpp      main.cpp
)

target_link_libraries(dmft ${CLUSTER_LIBRARIES} ${FOURIER_LIBRARIES} ${GREENS_LIBRARIES} ${CLUSTER_LIBRARIES} ${ALPS_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS dmft RUNTIME DESTINATION "bin/")
install(FILES cluster1d.xml cluster3d.xml cluster.xml DESTINATION "share/")
