add_library(alps_numeric_xml INTERFACE)
add_library(ALPS::numeric_xml ALIAS alps_numeric_xml)
set_target_properties(alps_numeric_xml PROPERTIES EXPORT_NAME numeric_xml)
target_link_libraries(alps_numeric_xml INTERFACE ALPS::numerics ALPS::xml)
target_sources(alps_numeric_xml INTERFACE
  FILE_SET numeric_xml_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/xml/matrix.hpp)

if(ALPS_BUILD_TESTING)
  add_executable(numeric_xml_contract "${PROJECT_SOURCE_DIR}/tests/cmake/consumer/numeric_xml.cpp")
  target_link_libraries(numeric_xml_contract PRIVATE ALPS::numeric_xml)
  alps_add_test(numeric_xml_contract)
  set_property(TEST numeric_xml_contract PROPERTY LABELS "numeric;xml;compatibility")
endif()
