# Module-owned implementation and public headers; linked into the common ALPS runtime.
target_sources(alps PRIVATE
  src/make_deprecated_parameters.cpp
  src/make_parameters_from_xml.cpp
  src/paramproxy.cpp
  src/params.cpp
  src/paramvalue.cpp)
target_sources(alps_headers INTERFACE
  FILE_SET params_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/detail/paramiterator.hpp
    include/alps/ngs/detail/paramproxy.hpp
    include/alps/ngs/detail/params_impl_base.hpp
    include/alps/ngs/detail/paramvalue.hpp
    include/alps/ngs/detail/paramvalue_reader.hpp
    include/alps/ngs/make_deprecated_parameters.hpp
    include/alps/ngs/make_parameters_from_xml.hpp
    include/alps/ngs/params.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
