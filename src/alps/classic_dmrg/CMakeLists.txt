# Classic (two-site, finite-system) DMRG on compiled models: Adrian Feiguin's
# dmtk engine (the Eigen-based rework of src/apps/dmrg/dmrg/dmtk, ported from
# the alpsx fork), driven by the numeric local terms of ALPS::models. The
# component depends on ALPS::models publicly and privately on Eigen (dense
# eigensolvers) and the platform BLAS (block products, as the legacy engine):
# every dmtk header is private to this library, so consumers of the installed
# SDK need neither Eigen nor BLAS headers. It does not use the XML catalog,
# legacy parameters, HDF5 or any other ALPS library.
add_library(alps_classic_dmrg src/classic.cpp)
target_link_libraries(alps_classic_dmrg
  PUBLIC ALPS::models
  PRIVATE Eigen3::Eigen ${ALPS_BLAS_TARGET})
target_include_directories(alps_classic_dmrg PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
alps_configure_component(alps_classic_dmrg classic_dmrg)
target_sources(alps_classic_dmrg PUBLIC
  FILE_SET classic_dmrg_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/dmrg/classic.hpp
    include/alps/dmrg/molecular.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
