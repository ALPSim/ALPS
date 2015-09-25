cmake_minimum_required(VERSION 3.0)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)
find_package(BlasClasses REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})

add_library (cluster
cluster_fourier_generator.cpp
CDMFT_transform.cpp          cluster_helper.cpp           
DCA_transform.cpp            
AFM_DCA_transform.cpp            
cluster_graph.cpp 
)

target_link_libraries(cluster ${ALPSCore_LIBRARIES})
