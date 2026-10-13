if(UNIX)
  file(RELATIVE_PATH ALPS_XML_FROM_BINDIR
    "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/alps/xml")
  configure_file(alps-xml.in alps-xml @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/alps-xml"
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
  # Retain the historical names with the same relocatable resource lookup.
  foreach(wrapper IN ITEMS convert2html convert2text extractgp extracthtml
      extracttext extractxmgr plot2gp plot2html plot2text plot2xmgr
      use_local_stylesheet)
    configure_file(${wrapper} ${wrapper} @ONLY NEWLINE_STYLE UNIX)
    install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${wrapper}"
      DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
  endforeach()
  install(PROGRAMS archivecat xslttransform
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
endif()

# txt2archive remains unbuilt/uninstalled.
