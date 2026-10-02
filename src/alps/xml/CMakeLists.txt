# Xml source ownership.
target_sources(alps PRIVATE
  src/parser/parser.C
  src/parser/xmlattributes.C
  src/parser/xmlhandler.C
  src/parser/xmlparser.C
  src/parser/xmlstream.C
  src/parser/xslt_path.C)

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
    include/alps/plot.h
    include/alps/xml.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME xml HEADER_SET xml_headers
  TARGETS alps
  DEPENDS configuration legacy_parameters utilities
)
