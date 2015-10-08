cmake_minimum_required(VERSION 3.0)

project (selfconsistency)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPSCore REQUIRED COMPONENTS hdf5 params)
find_package(GFTools REQUIRED)
#we also need the blas classes
find_package(BlasClasses REQUIRED)
find_package(Cluster REQUIRED)
find_package(Solver REQUIRED)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${SOLVER_INCLUDE_DIR})
include_directories(${GFTOOLS_INCLUDE_DIR})

add_library (selfconsistency
hilbert_transform.cpp        
selfconsistency.cpp          
#selfconsistency_CDMFT.cpp    selfconsistency_DCAPlus.cpp
selfconsistency_DCA.cpp
)

target_link_libraries(selfconsistency ${GFTOOLS_LIBRARIES} ${ALPSCore_LIBRARIES}
#${CLUSTER_LIBRARIES}
)
