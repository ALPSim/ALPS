include_directories(${SQLite_INCLUDE_DIR})

#
# parameter2xml
# 

add_executable(parameter2xml parameter2xml.C)
target_link_libraries(parameter2xml alps)
install(TARGETS parameter2xml RUNTIME DESTINATION bin)

#
# archive
#

if(SQLite_FOUND)
  add_executable(archive archive.cpp archive_index.cpp archive_node.cpp archive_plot.cpp
    archive_sqlite.cpp archive_xml.cpp)
  target_link_libraries(archive alps ${SQLite_LIBRARIES})
  install(TARGETS archive RUNTIME DESTINATION bin)
endif(SQLite_FOUND)

#
# lattice-preview and helper program
#

configure_file(config.py.in ${CMAKE_CURRENT_BINARY_DIR}/config.py)
if(WIN32 AND NOT UNIX)
  add_pi_executable(lattice-preview preview.py ${CMAKE_CURRENT_BINARY_DIR}/config.py license.py)
  file(GLOB pi_generated_files ${CMAKE_CURRENT_BINARY_DIR}/lattice-preview/*)
  install(FILES ${pi_generated_files} DESTINATION bin)
else(WIN32 AND NOT UNIX)
  configure_file(lattice-preview.in ${CMAKE_CURRENT_BINARY_DIR}/lattice-preview)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/lattice-preview DESTINATION bin)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.py preview.py license.py
    DESTINATION lib/python/alps)
endif(WIN32 AND NOT UNIX)
add_executable(lattice2xml lattice2xml.C)
target_link_libraries(lattice2xml alps)
install(TARGETS lattice2xml RUNTIME DESTINATION bin)
