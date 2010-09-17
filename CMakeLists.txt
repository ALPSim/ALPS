cmake_minimum_required(VERSION 2.8)
find_package(ALPS REQUIRED NO_SYSTEM_ENVIRONMENT_PATH)
include(${ALPS_USE_FILE})

project (DMFT)
include_directories("/Users/emanuel/alps_svn/")
add_executable(dmft 
2dsimpson.C clusterhelper.C fouriertransform.C selfconsistency.C
alps_solver.C clustertransformer.C hilberttransformer.C interpolator.C
auxiliaryfunctions.C externalsolver.C hirschfyesim.C main.C

../DMFT-Solvers/Werner/impurity.C  ../DMFT-Solvers/Werner/fourpoint.C

../DMFT-Solvers/Rubtsov/rubtsov.cpp
../DMFT-Solvers/Rubtsov/rubtsov_io.cpp
../DMFT-Solvers/Rubtsov/rubtsov_solver.cpp
../DMFT-Solvers/Rubtsov/rubtsov_auxiliary.cpp
../DMFT-Solvers/Rubtsov/rubtsov_splines.cpp
../DMFT-Solvers/Rubtsov/rubtsov_measurements.cpp
../DMFT-Solvers/Rubtsov/rubtsov_wang_landau.cpp 
../DMFT-Solvers/Rubtsov/rubtsov_fastupdate.cpp
../DMFT-Solvers/Rubtsov/rubtsov_model.cpp
../DMFT-Solvers/Rubtsov/rubtsov_worm.cpp
../DMFT-Solvers/Rubtsov/rubtsov_fourpoint.cpp
../DMFT-Solvers/Rubtsov/rubtsov_observables.cpp
../DMFT-Solvers/Rubtsov/rubtsov_selfenergy.cpp)
target_link_libraries(dmft ${ALPS_LIBRARIES})
