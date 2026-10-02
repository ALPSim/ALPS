# Existing command-line grammars are independent of simulation and input parsing.
add_library(alps_cli src/mcoptions.cpp src/parseargs.cpp)
add_library(ALPS::cli ALIAS alps_cli)
set_target_properties(alps_cli PROPERTIES
  EXPORT_NAME cli
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_cli PUBLIC alps_headers ALPS::utilities Boost::program_options)
alps_dependency_rpath(alps_cli)
alps_copy_runtime_dependencies(alps_cli)
generate_export_header(alps_cli EXPORT_MACRO_NAME ALPS_CLI_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/cli_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_CLI_STATIC_DEFINE)
endif()
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/cli_export.h")
target_sources(alps_headers INTERFACE
  FILE_SET cli_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/mcoptions.hpp include/alps/parseargs.hpp)

install(TARGETS alps_cli EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME cli HEADER_SET cli_headers
  TARGETS alps_cli
  DEPENDS utilities
  GENERATED_HEADERS alps/cli_export.h)
