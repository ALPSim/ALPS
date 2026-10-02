# Solvers source ownership.
target_sources(alps_headers INTERFACE
  FILE_SET solvers_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/solvers.hpp)

alps_register_module(NAME solvers HEADER_SET solvers_headers
  TARGETS alps
)
