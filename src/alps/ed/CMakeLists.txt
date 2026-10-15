# Exact diagonalization on compiled operators from ALPS::models: product-basis
# enumeration, sparse assembly, a dense LAPACK eigensolver and a deflated
# Lanczos driver over the in-tree IETL. The component depends on ALPS::models,
# Boost headers and LAPACK; it does not use the XML catalog, legacy parameters,
# HDF5 or any other ALPS library. IETL and the Boost numeric bindings are
# header-only build-time dependencies of the library sources.
add_library(alps_ed
  src/eigensolver.cpp
  src/lanczos.cpp
  src/product_basis.cpp
  src/product_problem.cpp
  src/sparse_from_local.cpp)
target_link_libraries(alps_ed
  PUBLIC ALPS::models
  PRIVATE ${ALPS_LAPACK_TARGET} ${ALPS_BLAS_TARGET})
target_include_directories(alps_ed PRIVATE
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/alps/ietl/include>")
alps_configure_component(alps_ed ed)
target_sources(alps_ed PUBLIC
  FILE_SET ed_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ed.hpp
    include/alps/ed/eigensolver.hpp
    include/alps/ed/lanczos.hpp
    include/alps/ed/product_basis.hpp
    include/alps/ed/product_problem.hpp
    include/alps/ed/sparse_from_local.hpp
    include/alps/ed/sparse_matrix.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
