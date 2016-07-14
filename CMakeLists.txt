cmake_minimum_required(VERSION 2.8)

project (DMFT)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(ALPSCore REQUIRED COMPONENTS hdf5 accumulators mc params gf)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -DNDEBUG -DBOOST_DISABLE_ASSERTS")
add_subdirectory(libcluster)
add_subdirectory(libselfconsistency)
add_subdirectory(tests)


find_package(LAPACK REQUIRED)
# FIXME: make this compiler-dependent
if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -ipo -no-prec-div -fp-model fast=2")
endif()

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})
include_directories("libcluster")
include_directories("libsolver")
include_directories("libselfconsistency")
include_directories("tests")

add_executable(dmft
main.cpp
)

enable_testing()

target_link_libraries(dmft selfconsistency cluster
${ALPSCore_LIBRARIES} ${LAPACK_LIBRARIES} 
)

install(TARGETS dmft RUNTIME DESTINATION "bin/")
FILE(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/blas_classes/*.h")
install(FILES ${files} DESTINATION "include/")


