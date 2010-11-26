cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
include_directories("/home/esg/alps_svn/")
add_executable(dmft 
2dsimpson.C clusterhelper.C fouriertransform.C selfconsistency.C
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
target_link_libraries(dmft ${ALPS_LIBRARIES} pthread)
