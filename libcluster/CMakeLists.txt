cmake_minimum_required(VERSION 3.0)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package(BlasClasses REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GFTools REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GFTOOLS_INCLUDE_DIR})

add_library (cluster
cluster_helper.cpp           
cluster_graph.cpp 
cluster_transform.cpp 
#DCA_transform.cpp            
#AFM_DCA_transform.cpp            
)

target_link_libraries(cluster  ${GFTOOLS_LIBRARIES} ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES})

include(EnableGtests)
UseGtest()
enable_testing()
include_directories(${GTEST_INCLUDE_DIR})
list(APPEND LINK_ALL cluster gtest_main)
add_gtest(cluster_helper_test)
add_gtest(cluster_transform_test)
