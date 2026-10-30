# Historical archive workflows remain available independently of alps-xml.
foreach(program IN ITEMS txt2archive xml2archive)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
set(archive_tools txt2archive xml2archive)

# Preserve the previous optional dependency; the SDK does not require SQLite.
find_package(SQLite3 QUIET)
if(SQLite3_FOUND)
  add_executable(archive archive.cpp archive_index.cpp archive_node.cpp
    archive_plot.cpp archive_sqlite.cpp archive_xml.cpp)
  target_link_libraries(archive PRIVATE ALPS::alps SQLite::SQLite3)
  list(APPEND archive_tools archive)
endif()
install(TARGETS ${archive_tools}
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
