cmake_minimum_required(VERSION 2.8)
project (AuxiliaryPrograms)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPSCore REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
#find_package(LAPACK REQUIRED)
include(find_blas_lapack)
find_package(BlasClasses REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
find_package(Cluster REQUIRED)
find_package(Selfconsistency REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -DNDEBUG -DBOOST_DISABLE_ASSERTS")


include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${SELFCONSISTENCY_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})

#add_executable(sigma_interpolate_sc sigma_interpolate_sc.cpp )
#
add_executable(sigma_interpolate_normal sigma_interpolate_normal.cpp )

target_link_libraries(sigma_interpolate_normal ${FOURIER_LIBRARIES} ${CLUSTER_LIBRARIES} ${SELFCONSISTENCY_LIBRARIES} ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES})
