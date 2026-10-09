# These entry points use the older Parameters and XML/parser implementation.
# Keep them with that runtime; the typed params library must not link back to it.
target_sources(alps PRIVATE
  params_from_file.cpp
  make_deprecated_parameters.cpp
  make_parameters_from_xml.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET params_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/params_from_file.hpp
        include/alps/ngs/make_deprecated_parameters.hpp
        include/alps/ngs/make_parameters_from_xml.hpp)
