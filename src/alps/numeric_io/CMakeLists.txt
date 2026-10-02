# Archive adapters depend on numerical types, never the reverse.
add_library(alps_numeric_io INTERFACE)
add_library(ALPS::numeric_io ALIAS alps_numeric_io)
set_target_properties(alps_numeric_io PROPERTIES EXPORT_NAME numeric_io)
target_link_libraries(alps_numeric_io INTERFACE ALPS::numerics ALPS::hdf5)
target_sources(alps_numeric_io INTERFACE
  FILE_SET numeric_io_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/hdf5/matrix.hpp
        include/alps/hdf5/numeric_vector.hpp)
alps_register_module(NAME numeric_io HEADER_SET numeric_io_headers
  HEADER_TARGETS alps_numeric_io
  DEPENDS numerics hdf5 utilities)
