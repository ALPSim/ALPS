cmake_minimum_required(VERSION 3.0)

project (green)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
#we also need the blas classes
find_package(BlasClasses REQUIRED)
find_package(Greens REQUIRED)
find_package(Fourier REQUIRED)

include(${ALPS_USE_FILE})
include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GREENS_INCLUDE_DIR})
include_directories(${FOURIER_INCLUDE_DIR})

add_library (cluster
cluster_fourier_generator.cpp
CDMFT_transform.cpp          cluster_helper.cpp           
DCA_transform.cpp            
AFM_transform.cpp            
alpsgraph.cpp 
)

target_link_libraries(cluster ${FOURIER_LIBRARIES} ${GREENS_LIBRARIES})
