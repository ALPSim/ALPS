# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPS::aleav2, ALEA v2: the ALEA library of the ALPS Core Libraries, copied from
# ALPSCore. ALPS's own ALEA, the default, is in ../alea.
# README.md records the source commit, the local changes and the license.
#
# Eigen is optional for ALPS as a whole, so ALEA v2 is built only when it is found.
# Check the version here: Eigen 5's version file rejects a "3.3" request, and
# Eigen 3's rejects a version range that reaches past 3.x.
find_package(Eigen3 CONFIG QUIET)
if(NOT Eigen3_FOUND OR Eigen3_VERSION VERSION_LESS 3.3)
  message(STATUS "Eigen 3.3 or newer not found: ALPS::aleav2 is not built")
  return()
endif()

add_library(alps_aleav2
  src/autocorr.cpp
  src/batch.cpp
  src/covariance.cpp
  src/fisher_f.cpp
  src/galois.cpp
  src/mean.cpp
  src/model.cpp
  src/propagation.cpp
  src/result.cpp
  src/testing.cpp
  src/transform.cpp
  src/util.cpp
  src/variance.cpp)
target_link_libraries(alps_aleav2 PUBLIC ALPS::configuration ALPS::hdf5 Eigen3::Eigen)
if(ALPS_ENABLE_MPI)
  # alps/alea/mpi.hpp reduces through ALPSCore's alps/utilities/mpi.hpp.
  target_link_libraries(alps_aleav2 PUBLIC ALPS::utilities MPI::MPI_CXX)
endif()
alps_configure_component(alps_aleav2 aleav2)
# ALPSCore's ALEA has no export macros; export everything as ALPSCore did.
set_target_properties(alps_aleav2 PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
target_sources(alps_aleav2 PUBLIC
  FILE_SET alea_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/alea.hpp
    include/alps/alea/autocorr.hpp
    include/alps/alea/batch.hpp
    include/alps/alea/bundle.hpp
    include/alps/alea/complex_op.hpp
    include/alps/alea/computed.hpp
    include/alps/alea/convert.hpp
    include/alps/alea/core.hpp
    include/alps/alea/covariance.hpp
    include/alps/alea/fisher_f.hpp
    include/alps/alea/hdf5.hpp
    include/alps/alea/internal/format.hpp
    include/alps/alea/internal/galois.hpp
    include/alps/alea/internal/joined.hpp
    include/alps/alea/internal/outer.hpp
    include/alps/alea/internal/pooling.hpp
    include/alps/alea/internal/util.hpp
    include/alps/alea/mean.hpp
    include/alps/alea/plugin/stream_serializer.hpp
    include/alps/alea/propagation.hpp
    include/alps/alea/result.hpp
    include/alps/alea/serialize.hpp
    include/alps/alea/testing.hpp
    include/alps/alea/transform.hpp
    include/alps/alea/transformer.hpp
    include/alps/alea/util.hpp
    include/alps/alea/util/model.hpp
    include/alps/alea/util/serializer.hpp
    include/alps/alea/var_strategy.hpp
    include/alps/alea/variance.hpp
    # ALPSCore keeps these in its utilities module; only ALEA v2 uses them here.
    include/alps/common/view.hpp
    include/alps/serialization/core.hpp
    include/alps/serialization/eigen.hpp
    include/alps/serialization/util.hpp)
if(ALPS_ENABLE_MPI)
  # It includes alps/utilities/mpi.hpp, which only MPI builds install.
  target_sources(alps_aleav2 PUBLIC FILE_SET alea_headers FILES include/alps/alea/mpi.hpp)
endif()

# A separate export, like the applications': consumers that do not ask for the
# aleav2 component do not need Eigen.
get_target_property(alea_header_sets alps_aleav2 INTERFACE_HEADER_SETS)
set(alea_header_install_args "")
foreach(header_set IN LISTS alea_header_sets)
  list(APPEND alea_header_install_args FILE_SET ${header_set}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers)
endforeach()
install(TARGETS alps_aleav2 EXPORT ALPSAleav2Targets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ${alea_header_install_args})
install(EXPORT ALPSAleav2Targets
  NAMESPACE ALPS::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps
  COMPONENT build)
install(FILES COPYRIGHT.TXT LICENSE.TXT
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/aleav2
  COMPONENT libraries)
set_property(GLOBAL PROPERTY ALPS_ALEAV2_BUILT ON)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
