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
2dsimpson.C clusterhelper.C fouriertransform.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
alps_solver.C AFM_transform.C DCA_transform.C hilberttransformer.C interpolator.C
auxiliaryfunctions.C externalsolver.C main.C)

add_executable(kspace_interpolate interpolate.C clusterhelper.C
fouriertransform.C auxiliaryfunctions.C
AFM_transform.C DCA_transform.C)
target_link_libraries(kspace_interpolate ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
add_executable(sigma_interpolate_sc sigma_interpolate_sc.C
auxiliaryfunctions.C AFM_transform.C DCA_transform.C clusterhelper.C
fouriertransform.C)

add_executable(sigma_interpolate sigma_interpolate_normal.C
auxiliaryfunctions.C AFM_transform.C DCA_transform.C clusterhelper.C
fouriertransform.C)

add_executable(bg_parameter bg_parameter.C)

add_executable(hdf5_changeval hdf5_changeval.C)

add_executable(AFM_test AFM_test.C AFM_transform.C DCA_transform.C clusterhelper.C auxiliaryfunctions.C fouriertransform.C)

add_executable(AFM_4site AFM_4site.C AFM_transform.C DCA_transform.C clusterhelper.C auxiliaryfunctions.C fouriertransform.C)


target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
#target_link_libraries(energy ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
target_link_libraries(sigma_interpolate gsl gslcblas ${ALPS_LIBRARIES}
${LAPACK_LIBRARY} pthread)
target_link_libraries(sigma_interpolate_sc gsl gslcblas ${ALPS_LIBRARIES}
${LAPACK_LIBRARY} pthread )
target_link_libraries(kspace_interpolate ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})

target_link_libraries(bg_parameter ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
target_link_libraries(hdf5_changeval ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})

target_link_libraries(AFM_test ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})

target_link_libraries(AFM_4site ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
