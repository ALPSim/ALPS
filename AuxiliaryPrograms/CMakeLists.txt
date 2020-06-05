cmake_minimum_required(VERSION 2.8)
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()
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

foreach(f normal sc)
    add_executable(sigma_interpolate_${f} sigma_interpolate_${f}.cpp)
    target_link_libraries(sigma_interpolate_${f} ${FOURIER_LIBRARIES} ${CLUSTER_LIBRARIES} ${SELFCONSISTENCY_LIBRARIES} ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} ${GSL_LIBRARIES} ${MPI_CXX_LIBRARIES})
endforeach()