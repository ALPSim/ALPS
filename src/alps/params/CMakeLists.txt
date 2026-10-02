# Typed parameters have no link dependency on the legacy parameter/parser runtime.
add_library(alps_params
  src/paramproxy.cpp
  src/params.cpp
  src/paramvalue.cpp)
add_library(ALPS::params ALIAS alps_params)
set_target_properties(alps_params PROPERTIES
  EXPORT_NAME params
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_params PUBLIC alps_headers ALPS::hdf5 Boost::serialization)
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps_params PUBLIC MPI::MPI_CXX Boost::mpi)
endif()
alps_dependency_rpath(alps_params)
alps_copy_runtime_dependencies(alps_params)
generate_export_header(alps_params EXPORT_MACRO_NAME ALPS_PARAMS_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/params_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_PARAMS_STATIC_DEFINE)
endif()
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/params_export.h")
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

install(TARGETS alps_params EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory(adapters)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
