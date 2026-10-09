if(UNIX)
  file(RELATIVE_PATH ALPS_XML_FROM_BINDIR
    "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/alps/xml")
  configure_file(alps-xml.in alps-xml @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/alps-xml"
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
endif()

# txt2archive and the historical shell wrappers remain unbuilt/uninstalled.

# Retain the historical command names alongside the new frontend.
if(UNIX)
  install(PROGRAMS archivecat convert2html convert2text extractgp extracthtml
    extracttext extractxmgr plot2gp plot2html plot2text plot2xmgr
    use_local_stylesheet xslttransform
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
endif()
