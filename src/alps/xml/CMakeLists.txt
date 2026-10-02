# XML parsing and output are independent of the simulation runtime.
add_library(alps_xml
  src/parser/parser.C
  src/parser/xmlattributes.C
  src/parser/xmlhandler.C
  src/parser/xmlparser.C
  src/parser/xmlstream.C
  src/parser/xslt_path.C)
add_library(ALPS::xml ALIAS alps_xml)
set_target_properties(alps_xml PROPERTIES
  EXPORT_NAME xml
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_xml PUBLIC alps_headers Boost::filesystem
  PRIVATE Boost::regex)
alps_dependency_rpath(alps_xml)
alps_copy_runtime_dependencies(alps_xml)
generate_export_header(alps_xml EXPORT_MACRO_NAME ALPS_XML_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/xml_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_XML_STATIC_DEFINE)
endif()
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/xml_export.h")

target_sources(alps_headers INTERFACE
  FILE_SET xml_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/parser/parser.h
    include/alps/parser/xmlattributes.h
    include/alps/parser/xmlhandler.h
    include/alps/parser/xmlparser.h
    include/alps/parser/xmlstream.h
    include/alps/parser/xslt_path.h
    include/alps/xml.h)

install(TARGETS alps_xml EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME xml HEADER_SET xml_headers
  TARGETS alps_xml
  DEPENDS configuration utilities
  GENERATED_HEADERS alps/xml_export.h
)
