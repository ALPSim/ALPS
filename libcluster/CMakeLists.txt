#cmake_minimum_required(VERSION 3.0)
cmake_minimum_required(VERSION 2.8.12)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
# find_package(MPI REQUIRED) # already in ALPSCore
find_package(BlasClasses REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GFTools REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GFTools_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MPI_CXX_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

add_library (cluster
cluster_helper.cpp           
cluster_graph.cpp 
cluster_transform.cpp 
)

target_link_libraries(cluster  ${GFTools_LIBRARIES} ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} #${MPI_CXX_LIBRARIES}
)

include(EnableGtests)
UseGtest()
enable_testing()
include_directories(${GTEST_INCLUDE_DIR})
list(APPEND LINK_ALL cluster gtest_main)
add_gtest(cluster_helper_test)
add_gtest(cluster_transform_test)
