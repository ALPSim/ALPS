include_directories(${SQLite_INCLUDE_DIR})

add_executable(parameter2xml parameter2xml.C)
target_link_libraries(parameter2xml alps)
install(TARGETS parameter2xml RUNTIME DESTINATION bin)

if(SQLite_FOUND)
  add_executable(archive archive.cpp archive_index.cpp archive_node.cpp archive_plot.cpp
    archive_sqlite.cpp archive_xml.cpp)
  target_link_libraries(archive alps ${SQLite_LIBRARIES})
  install(TARGETS archive RUNTIME DESTINATION bin)
endif(SQLite_FOUND)
