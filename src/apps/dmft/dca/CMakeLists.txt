# DCA: dynamical cluster approximation self-consistency loop, formerly
# CQMP/DMFT. The impurity problem is solved by an external solver executable
# that exchanges HDF5 files with the loop.
# Copyright (C) 2005-2026 by Emanuel Gull
# SPDX-License-Identifier: MIT

# The cluster transformers store Green's functions in ALPS::gf, which needs Eigen.
if(NOT TARGET ALPS::gf)
  message(STATUS "ALPS::gf is not built (Eigen not found): dca is not built")
  return()
endif()

# The cluster geometry, k-space transforms and Dyson equations, the parameters
# of the self-consistency loop and the loop itself, shared by the dca program,
# the self-energy interpolation tools and the tests.
add_library(alps_dca STATIC
    libcluster/cluster_graph.cpp
    libcluster/cluster_helper.cpp
    libcluster/cluster_transform.cpp
    libselfconsistency/selfconsistency.cpp
    libselfconsistency/selfconsistency_DCA.cpp)
target_include_directories(alps_dca PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/libcluster"
    "${CMAKE_CURRENT_SOURCE_DIR}/libselfconsistency"
    "${CMAKE_CURRENT_SOURCE_DIR}/libsolver"
    "${CMAKE_CURRENT_SOURCE_DIR}/blas_classes")
target_link_libraries(alps_dca PUBLIC ALPS::alps ALPS::gf)
# alps/gf/flagcheck.hpp asks consumers to choose its debug checks, as ALPS::gf does.
target_compile_definitions(alps_dca PUBLIC
  "$<IF:$<CONFIG:Debug>,ALPS_GF_DEBUG,BOOST_DISABLE_ASSERTS>")

add_executable(dca main.cpp)
target_link_libraries(dca PRIVATE alps_dca)
alps_install_applications(dca)

# Cluster definitions read through dca.LATTICE_LIBRARY.
install(FILES libcluster/cluster.xml libcluster/cluster1d.xml libcluster/cluster3d.xml
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/dca COMPONENT applications)

# Tools that interpolate a converged self-energy onto another cluster, to start
# a new calculation from it. They need GSL for the splines.
find_package(GSL QUIET)
if(GSL_FOUND)
  foreach(kind IN ITEMS normal sc)
    add_executable(sigma_interpolate_${kind} AuxiliaryPrograms/sigma_interpolate_${kind}.cpp)
    target_link_libraries(sigma_interpolate_${kind} PRIVATE alps_dca GSL::gsl)
  endforeach()
  alps_install_applications(sigma_interpolate_normal sigma_interpolate_sc)
else()
  message(STATUS "GSL not found: sigma_interpolate_normal and sigma_interpolate_sc are not built")
endif()

if(ALPS_BUILD_TESTING)
  foreach(test IN ITEMS cluster_helper_test cluster_transform_test
                        double_cell_transform_test normal_state_transform_test)
    add_executable(dca_${test} tests/${test}.cpp)
    target_link_libraries(dca_${test} PRIVATE alps_dca)
    target_compile_definitions(dca_${test} PRIVATE
      DCA_CLUSTER_LIBRARY="${CMAKE_CURRENT_SOURCE_DIR}/libcluster/cluster.xml")
    alps_add_gtest(dca_${test} COMPONENT dmft)
  endforeach()
endif()
