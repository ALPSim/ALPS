# Graph source ownership.
target_sources(alps_headers INTERFACE
  FILE_SET graph_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/graph/canonical_graph.hpp
    include/alps/graph/canonical_properties.hpp
    include/alps/graph/canonical_properties_traits.hpp
    include/alps/graph/detail/assert_helpers.hpp
    include/alps/graph/detail/canonical_properties_impl.hpp
    include/alps/graph/detail/helper_functions.hpp
    include/alps/graph/detail/lattice_constant_impl.hpp
    include/alps/graph/detail/shared_queue.hpp
    include/alps/graph/is_embeddable.hpp
    include/alps/graph/lattice_constant.hpp
    include/alps/graph/lattice_constant_debug.hpp
    include/alps/graph/subgraph_generator.hpp
    include/alps/graph/subgraphs.hpp
    include/alps/graph/utils.hpp
    include/alps/graph/vertices_of_cell.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME graph HEADER_SET graph_headers
  TARGETS alps
  DEPENDS lattice numerics utilities
  DEBT_FILE architecture-debt.json
)
