# These entry points use the older Parameters and XML/parser implementation.
# Keep them with that runtime; the typed params library must not link back to it.
target_sources(alps PRIVATE
  params_from_file.cpp
  make_deprecated_parameters.cpp
  make_parameters_from_xml.cpp)
