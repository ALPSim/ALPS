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
target_link_libraries(alps_xml PUBLIC ALPS::configuration "$<COMPILE_ONLY:ALPS::utilities>" Boost::filesystem
  PRIVATE Boost::regex)
alps_dependency_rpath(alps_xml)
alps_copy_runtime_dependencies(alps_xml)
generate_export_header(alps_xml EXPORT_MACRO_NAME ALPS_XML_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/xml_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_xml PUBLIC ALPS_XML_STATIC_DEFINE)
endif()
target_sources(alps_xml PUBLIC
  FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include"
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/xml_export.h")

target_sources(alps_xml PUBLIC
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

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME xml HEADER_SET xml_headers
  HEADER_TARGETS alps_xml
  TARGETS alps_xml
  DEPENDS configuration utilities
  GENERATED_HEADERS alps/xml_export.h
)
