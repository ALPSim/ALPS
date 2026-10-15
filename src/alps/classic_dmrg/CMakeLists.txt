# Classic (two-site, finite-system) DMRG on compiled models: Adrian Feiguin's
# dmtk engine (the Eigen-based rework of src/apps/dmrg/dmrg/dmtk, ported from
# the alpsx fork), driven by the numeric local terms of ALPS::models. The
# component depends on ALPS::models publicly and on Eigen privately: every
# dmtk header is private to this library, so consumers of the installed SDK
# need neither Eigen nor LAPACK. It does not use the XML catalog, legacy
# parameters, HDF5 or any other ALPS library.
add_library(alps_classic_dmrg src/classic.cpp)
target_link_libraries(alps_classic_dmrg
  PUBLIC ALPS::models
  PRIVATE Eigen3::Eigen)
target_include_directories(alps_classic_dmrg PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
alps_configure_component(alps_classic_dmrg classic_dmrg)
target_sources(alps_classic_dmrg PUBLIC
  FILE_SET classic_dmrg_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/dmrg/classic.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
