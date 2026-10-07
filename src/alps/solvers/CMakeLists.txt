add_library(alps_solver_headers INTERFACE)
add_library(ALPS::solver_headers ALIAS alps_solver_headers)
set_target_properties(alps_solver_headers PROPERTIES EXPORT_NAME solver_headers)
target_link_libraries(alps_solver_headers INTERFACE ALPS::configuration)
# Solvers source ownership.
target_sources(alps_solver_headers INTERFACE
  FILE_SET solvers_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/solvers.hpp)
