cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
#include_directories("/home/esg/alps_svn/")
include_directories("/Users/emanuel/alps_svn/")
include_directories("/Users/egull/alps_svn/")
include_directories("/home/gull/alps_svn/")
include_directories("/home2/egull/alps_svn")
#include_directories("/opt/local/include")
include_directories("/global/homes/g/gull/alps_svn/")
include_directories("/home2/egull/software/gsl-1.9/include")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ffast-math -march=native -O3 -DUSE_MPI -DNDEBUG")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -xHOST -O3 -DNDEBUG -mkl=sequential")

link_directories("/opt/local/lib")
link_directories("/home2/egull/software/gsl-1.9/lib")


add_executable(dmft 
clusterhelper.C fouriertransform.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
AFM_transform.C DCA_transform.C hilberttransformer.C 
auxiliaryfunctions.C externalsolver.C main.C)

target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
