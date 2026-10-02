if(UNIX)
  file(RELATIVE_PATH ALPS_XML_FROM_BINDIR
    "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/alps/xml")
  configure_file(alps-xml.in alps-xml @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/alps-xml"
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
endif()

# txt2archive and the historical shell wrappers remain unbuilt/uninstalled.
alps_register_module(NAME tools_xml
  EXTRA_FILES txt2archive.C
  DEPENDS xml)
