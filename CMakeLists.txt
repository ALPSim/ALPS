cmake_minimum_required(VERSION 2.8)

project (DMFT)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params)

add_subdirectory(libcluster)
add_subdirectory(libselfconsistency)

find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(GFTools REQUIRED)

include(UseGtest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHost -fast -march=native -ffast-math -DNDEBUG -DBOOST_DISABLE_ASSERTS")

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${GFTOOLS_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories("libcluster")
include_directories("libsolver")
include_directories("libselfconsistency")

add_executable(dmft
main.cpp
)

target_link_libraries(dmft selfconsistency cluster
${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} ${GFTOOLS_LIBRARIES} ${MPI_CXX_LIBRARIES})

install(TARGETS dmft RUNTIME DESTINATION "bin/")
install(FILES cluster1d.xml cluster3d.xml cluster.xml DESTINATION "share/")
