# XML parsing and output are independent of the simulation runtime.
add_library(alps_xml
  src/parser/parser.C
  src/parser/xmlattributes.C
  src/parser/xmlhandler.C
  src/parser/xmlparser.C
  src/parser/xmlstream.C
  src/parser/xslt_path.C)
target_link_libraries(alps_xml PUBLIC ALPS::configuration "$<COMPILE_ONLY:ALPS::utilities>" Boost::filesystem
  PRIVATE Boost::regex)
alps_configure_component(alps_xml xml)

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
