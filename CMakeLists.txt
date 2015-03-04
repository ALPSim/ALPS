cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(UseGtest)

project (DMFT)
include_directories("$ENV{HOME}/alps_svn/")
include_directories("/usr/include/gsl")
include_directories("${CMAKE_SOURCE_DIR}/gtest/include")
#SET(CMAKE_CXX_COMPILER "CC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -march=native -O3")

link_directories("/opt/local/lib")
link_directories("/usr/lib")
link_directories("${CMAKE_SOURCE_DIR}/gtest")

add_executable(dmft 
clusterhelper.C fouriertransform.C fourier_generator.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
AFM_transform.C DCA_transform.C CDMFT_transform.C hilberttransformer.C 
auxiliaryfunctions.C externalsolver.C main.C)

add_executable(runTests test.cpp fouriertransform.C)
target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
target_link_libraries(runTests ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread gtest gtest_main)
