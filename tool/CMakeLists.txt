include_directories(${SQLite_INCLUDE_DIR})

#
# basic C++ tools
# 

if (ALPS_FOR_VISTRAILS)
  if (WIN32 AND NOT UNIX)
    install(PROGRAMS msxsl.exe DESTINATION ${VISTRAILS_APP_NAME})
  endif(WIN32 AND NOT UNIX)
else(ALPS_FOR_VISTRAILS)
  set (TOOLS parameter2xml convert2xml compactrun lattice2xml pevaluate pconfig poutput printgraph txt2archive xml2archive)
  foreach(program ${TOOLS})
    add_executable(${program} ${program}.C)
    target_link_libraries(${program} alps)
    install(TARGETS ${program} RUNTIME DESTINATION bin) 
  endforeach(program)

  if(SQLite_FOUND)
    add_executable(archive archive.cpp archive_index.cpp archive_node.cpp archive_plot.cpp
      archive_sqlite.cpp archive_xml.cpp)
    target_link_libraries(archive alps ${SQLite_LIBRARIES})
    install(TARGETS archive RUNTIME DESTINATION bin)
  endif(SQLite_FOUND)

#
# legacy shell scripts that stay around for backwards compatibility until the new solutions are there
#

  if(UNIX AND NOT WIN32)
    set (SCRIPTS_TO_CONFIGURE convert2html convert2text extractgp extracthtml extractmpl extracttext extractxmgr plot2gp plot2html plot2mpl plot2text plot2xmgr)
    foreach(script ${SCRIPTS_TO_CONFIGURE})
      configure_file(${script}.in ${CMAKE_CURRENT_BINARY_DIR}/${script})
      install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin)
    endforeach(script)
    install(PROGRAMS archivecat xslttransform DESTINATION bin)
  else(UNIX AND NOT WIN32)
    install(PROGRAMS msxsl.exe DESTINATION bin)
  endif(UNIX AND NOT WIN32)

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
endif(ALPS_FOR_VISTRAILS)

#
# Analytical continuation with MaxEnt
#

if(LAPACK_FOUND AND NOT ALPS_FOR_VISTRAILS)
  add_definitions(${LAPACK_DEFINITIONS} -DBIND_FORTRAN_LOWERCASE_UNDERSCORE)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(maxent maxent.cpp maxent_helper.cpp maxent_simulation.cpp maxent_parms.cpp)
  target_link_libraries(maxent alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  install(TARGETS maxent RUNTIME DESTINATION bin)
endif(LAPACK_FOUND AND NOT ALPS_FOR_VISTRAILS)
