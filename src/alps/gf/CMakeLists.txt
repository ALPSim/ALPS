# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# ALPS::gf: the Green's function library of the ALPS Core Libraries, copied from
# ALPSCore with the tensors and type traits it needs. README.md records the
# source commit and the local changes.
if(NOT ALPS_HAVE_EIGEN)
  return()
endif()

add_library(alps_gf src/mesh.cpp)
target_link_libraries(alps_gf PUBLIC ALPS::configuration ALPS::hdf5 Eigen3::Eigen)
if(ALPS_ENABLE_MPI)
  # gf broadcasts through ALPSCore's alps/utilities/mpi.hpp.
  target_link_libraries(alps_gf PUBLIC ALPS::utilities MPI::MPI_CXX)
endif()
alps_configure_component(alps_gf gf)
# ALPSCore's gf has no export macros; export everything as ALPSCore did.
set_target_properties(alps_gf PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
# ALPSCore's build defined these in its own compiler flags, without passing them
# on to consumers; alps/gf/flagcheck.hpp warns when neither is defined.
set(ALPS_GF_ASSERT_DEFINITION "$<IF:$<CONFIG:Debug>,ALPS_GF_DEBUG,BOOST_DISABLE_ASSERTS>")
target_compile_definitions(alps_gf PRIVATE ${ALPS_GF_ASSERT_DEFINITION})
target_sources(alps_gf PUBLIC
  FILE_SET gf_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/gf/flagcheck.hpp
    include/alps/gf/fourier.hpp
    include/alps/gf/gf.hpp
    include/alps/gf/gf_base.hpp
    include/alps/gf/grid.hpp
    include/alps/gf/mesh.hpp
    include/alps/gf/mesh/chebyshev_mesh.hpp
    include/alps/gf/mesh/frequency_meshes.hpp
    include/alps/gf/mesh/index.hpp
    include/alps/gf/mesh/index_meshes.hpp
    include/alps/gf/mesh/legendre_mesh.hpp
    include/alps/gf/mesh/mesh_base.hpp
    include/alps/gf/mesh/numerical_mesh.hpp
    include/alps/gf/mesh/time_meshes.hpp
    include/alps/gf/piecewise_polynomial.hpp
    include/alps/gf/tail.hpp
    include/alps/gf/tail_base.hpp
    # ALPSCore keeps these in its utilities and hdf5 modules; only gf uses them here.
    include/alps/hdf5/tensor.hpp
    include/alps/numeric/tensors.hpp
    include/alps/numeric/tensors/data_storage.hpp
    include/alps/numeric/tensors/data_view.hpp
    include/alps/numeric/tensors/tensor_base.hpp
    include/alps/numeric/tensors/view.hpp
    include/alps/type_traits/are_all_integrals.hpp
    include/alps/type_traits/common.hpp
    include/alps/type_traits/index_sequence.hpp
    include/alps/type_traits/tuple_traits.hpp)
if(ALPS_ENABLE_MPI)
  # It includes alps/utilities/mpi.hpp, which only MPI builds install, and gf
  # includes it only when ALPS_HAVE_MPI is defined.
  target_sources(alps_gf PUBLIC FILE_SET gf_headers FILES include/alps/gf/mpi_bcast.hpp)
endif()
alps_install_alpscore_component(alps_gf gf)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
