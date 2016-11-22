
add_executable(unpack_measurement_terms unpack_measurement_terms.cpp)
target_link_libraries(unpack_measurement_terms ${MYAPP_LIBRARIES})
add_alps_test(unpack_measurement_terms)
set_tests_properties(unpack_measurement_terms
  PROPERTIES ENVIRONMENT "ALPS_XML_PATH=${CMAKE_BINARY_DIR}/lib/xml"
)

add_executable(product_operator product_operator.cpp)
target_link_libraries(product_operator ${MYAPP_LIBRARIES})
add_alps_test(product_operator)
set_tests_properties(product_operator
  PROPERTIES ENVIRONMENT "ALPS_XML_PATH=${CMAKE_BINARY_DIR}/lib/xml"
)

add_executable(parse_measurements parse_measurements.cpp)
target_link_libraries(parse_measurements ${MYAPP_LIBRARIES})
add_alps_test(parse_measurements)

add_executable(test_measurements measurements.cpp)
target_link_libraries(test_measurements ${MYAPP_LIBRARIES})
add_test(NAME test_mps_measurements COMMAND test_measurements "--"
 ${CMAKE_CURRENT_SOURCE_DIR}/custom_lattice.xml
 ${CMAKE_CURRENT_SOURCE_DIR}/measurements_ref.h5
)
set_tests_properties(test_mps_measurements
  PROPERTIES ENVIRONMENT "ALPS_XML_PATH=${CMAKE_BINARY_DIR}/lib/xml"
)
