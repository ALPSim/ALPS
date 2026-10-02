# Legacy parameters source ownership.
target_sources(alps PRIVATE
  src/parameter/parameter.C
  src/parameter/parameterlist.C
  src/parameter/parameters.C)

target_sources(alps_headers INTERFACE
  FILE_SET legacy_parameters_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/parameter.h
    include/alps/parameter/parameter.h
    include/alps/parameter/parameter_p.h
    include/alps/parameter/parameterlist.h
    include/alps/parameter/parameterlist_p.h
    include/alps/parameter/parameters.h
    include/alps/parameter/parameters_p.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME legacy_parameters HEADER_SET legacy_parameters_headers
  TARGETS alps
  DEPENDS configuration expression hdf5 osiris utilities xml
)
