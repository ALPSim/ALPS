cmake_minimum_required(VERSION 3.1)

project (selfconsistency)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../cmake)

# we need old ALPS
find_package(ALPSCore REQUIRED COMPONENTS hdf5 params gf)
find_package(MPI REQUIRED)
#we also need the blas classes
find_package(BlasClasses REQUIRED)
if (NOT TARGET cluster)
  find_package(Cluster REQUIRED)
else()
  set(CLUSTER_LIBRARIES "cluster")
  # FIXME: use INTERFACE_INCLUDE_DIRECTORIES in ../libcluster/CMakeLists.txt instead
  get_target_property(CLUSTER_INCLUDE_DIR cluster INCLUDE_DIRECTORIES)
  # message(STATUS "CLUSTER_INCLUDE_DIR=${CLUSTER_INCLUDE_DIR}")
endif()

if (NOT TARGET solver)
find_package(Solver REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DBOOST_DISABLE_ASSERTS")

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${SOLVER_INCLUDE_DIR})

add_library (selfconsistency
hilbert_transform.cpp        
selfconsistency.cpp          
#selfconsistency_CDMFT.cpp    selfconsistency_DCAPlus.cpp
selfconsistency_DCA.cpp
)

target_link_libraries(selfconsistency ${ALPSCore_LIBRARIES} ${CLUSTER_LIBRARIES}
 ${MPI_CXX_LIBRARIES}
)

install(TARGETS selfconsistency
   ARCHIVE DESTINATION "lib/"
   LIBRARY DESTINATION "lib/"
   COMPONENT library)

FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${files} DESTINATION include)
