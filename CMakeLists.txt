cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
include_directories("$ENV{HOME}/alps_svn/")
#include_directories("$ENV{GSL_DIR}/include")
#SET(CMAKE_CXX_COMPILER "CC")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -DNDEBUG -framework accelerate")

#link_directories("/opt/local/lib")
link_directories("/home/leblanc/software/gsl-1.14/install/lib/libgsl.so")



add_executable(dmft 
clusterhelper.C fouriertransform.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
AFM_transform.C DCA_transform.C CDMFT_transform.C hilberttransformer.C 
auxiliaryfunctions.C externalsolver.C main.C)

target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
