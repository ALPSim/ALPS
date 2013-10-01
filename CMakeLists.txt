cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
include_directories("$ENV{HOME}/alps_svn/")
include_directories("$ENV{GSL_DIR}/include")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -ipo -xHOST -openmp -O3 -DUSE_MPI -DNDEBUG -DMPICH_IGNORE_CXX_SEEK")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DMPICH_IGNORE_CXX_SEEK")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DUSE_MPI -DNDEBUG -DMPICH_IGNORE_CXX_SEEK")

link_directories("/opt/local/lib")
link_directories("/usr/lib/alps/")
link_directories("$ENV{GSL_DIR}/lib")

#link_libraries("/usr/lib/alps/libalpslli.so.0")
#link_libraries("/usr/lib/alps/libalpsutil.so.0")


add_executable(dmft 
2dsimpson.C clusterhelper.C fouriertransform.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
alps_solver.C clustertransformer.C hilberttransformer.C interpolator.C
auxiliaryfunctions.C externalsolver.C main.C)

#add_executable(sigma_interpolate_sc sigma_interpolate_sc.C
#auxiliaryfunctions.C clustertransformer.C clusterhelper.C
#fouriertransform.C)

#add_executable(sigma_interpolate sigma_interpolate_normal.C
#auxiliaryfunctions.C clustertransformer.C clusterhelper.C
#fouriertransform.C)

#add_executable(bg_parameter bg_parameter.C)

#add_executable(hdf5_changeval hdf5_changeval.C)


target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
#target_link_libraries(sigma_interpolate ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} -lgsl -lgslcblas )
#target_link_libraries(sigma_interpolate_sc ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} -lgsl -lgslcblas)

#target_link_libraries(bg_parameter ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
#target_link_libraries(hdf5_changeval ${ALPS_LIBRARIES} ${LAPACK_LIBRARY})
