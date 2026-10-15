# Programmatic lattice and model construction. The component depends only on
# the configuration interface (Boost headers); it does not use the XML catalog,
# legacy parameters, HDF5 or any other ALPS library.
add_library(alps_models
  src/basis.cpp
  src/builtin.cpp
  src/finite_graph.cpp
  src/geometry.cpp
  src/hamiltonian.cpp
  src/lattice.cpp
  src/physical_basis.cpp
  src/site.cpp
  src/terms.cpp)
target_link_libraries(alps_models PUBLIC ALPS::configuration)
alps_configure_component(alps_models models)
target_sources(alps_models PUBLIC
  FILE_SET models_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/models.hpp
    include/alps/models/builtin.hpp
    include/alps/models/compiled.hpp
    include/alps/models/finite_graph.hpp
    include/alps/models/geometry.hpp
    include/alps/models/hamiltonian.hpp
    include/alps/models/lattice.hpp
    include/alps/models/physical_basis.hpp
    include/alps/models/site.hpp
    include/alps/models/terms.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
