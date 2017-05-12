#cmake_minimum_required(VERSION 3.0)
cmake_minimum_required(VERSION 2.8.12)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params gf)
find_package(BlasClasses REQUIRED)
include(find_blas_lapack)
find_package(MPI REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DBOOST_DISABLE_ASSERTS")

add_library (cluster
cluster_helper.cpp           
cluster_graph.cpp 
cluster_transform.cpp 
)

target_link_libraries(cluster ${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} ${MPI_CXX_LIBRARIES}
)

install(TARGETS cluster 
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} DESTINATION include)
install(FILES cluster1d.xml cluster3d.xml cluster.xml DESTINATION "share/")

#include(EnableGtests)
#UseGtest()
#enable_testing()
#include_directories(${GTEST_INCLUDE_DIR})
#list(APPEND LINK_ALL cluster gtest_main)
#add_gtest(cluster_helper_test)
#add_gtest(cluster_transform_test)
