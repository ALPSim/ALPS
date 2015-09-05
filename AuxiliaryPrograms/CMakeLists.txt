cmake_minimum_required(VERSION 2.8)
project (AuxiliaryPrograms)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
find_package(LAPACK REQUIRED)
find_package(BlasClasses REQUIRED)
find_package(Greens REQUIRED)
find_package(Fourier REQUIRED)
find_package(Cluster REQUIRED)
find_package(GSL REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GREENS_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})

add_executable(sigma_interpolate_sc sigma_interpolate_sc.cpp )

add_executable(sigma_interpolate_normal sigma_interpolate_normal.cpp )

target_link_libraries(sigma_interpolate_normal ${GREENS_LIBRARIES} ${FOURIER_LIBRARIES} ${CLUSTER_LIBRARIES} ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} ${GSL_LIBRARIES})
target_link_libraries(sigma_interpolate_sc  ${GREENS_LIBRARIES} ${FOURIER_LIBRARIES} ${CLUSTER_LIBRARIES} ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} ${GSL_LIBRARIES})
