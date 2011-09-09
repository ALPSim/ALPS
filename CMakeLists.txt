cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
#include_directories("/home/esg/alps_svn/")
include_directories("/Users/emanuel/alps_svn/")
#include_directories("/opt/local/include")
include_directories("/cluster/home/phys/gullc/alps_svn")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=core2 -O3 -DUSE_MPI -DNDEBUG")

link_directories("/opt/local/lib")


add_executable(dmft 
2dsimpson.C clusterhelper.C fouriertransform.C selfconsistency.C
selfconsistency_DCA.C selfconsistency_CDMFT.C
alps_solver.C clustertransformer.C hilberttransformer.C interpolator.C
auxiliaryfunctions.C externalsolver.C hirschfyesim.C main.C

../DMFT-solvers/Werner/impurity.C  ../DMFT-solvers/Werner/fourpoint.C

../DMFT-solvers/Rubtsov/rubtsov.cpp
../DMFT-solvers/Rubtsov/rubtsov_io.cpp
../DMFT-solvers/Rubtsov/rubtsov_solver.cpp
../DMFT-solvers/Rubtsov/rubtsov_auxiliary.cpp
../DMFT-solvers/Rubtsov/rubtsov_splines.cpp
../DMFT-solvers/Rubtsov/rubtsov_measurements.cpp
../DMFT-solvers/Rubtsov/rubtsov_wang_landau.cpp 
../DMFT-solvers/Rubtsov/rubtsov_fastupdate.cpp
../DMFT-solvers/Rubtsov/rubtsov_model.cpp
../DMFT-solvers/Rubtsov/rubtsov_worm.cpp
../DMFT-solvers/Rubtsov/rubtsov_fourpoint.cpp
../DMFT-solvers/Rubtsov/rubtsov_observables.cpp
../DMFT-solvers/Rubtsov/rubtsov_selfenergy.cpp)

add_executable(energy ../DMFT-solvers/Rubtsov/energy.cpp
../DMFT-Solvers/Rubtsov/obs_auxiliary.cpp clustertransformer.C
hilberttransformer.C clusterhelper.C fouriertransform.C 2dsimpson.C
auxiliaryfunctions.C
)

add_executable(kspace_interpolate interpolate.C clusterhelper.C
fouriertransform.C auxiliaryfunctions.C
clustertransformer.C)
target_link_libraries(kspace_interpolate ${ALPS_LIBRARIES})
add_executable(sigma_interpolate sigma_interpolate.C
auxiliaryfunctions.C clustertransformer.C clusterhelper.C
fouriertransform.C)
target_link_libraries(dmft ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
target_link_libraries(energy ${ALPS_LIBRARIES} ${LAPACK_LIBRARY} pthread)
target_link_libraries(sigma_interpolate ${ALPS_LIBRARIES}
${LAPACK_LIBRARY} gsl pthread)

add_executable(bg_parameter bg_parameter.C)
target_link_libraries(bg_parameter ${ALPS_LIBRARIES})

