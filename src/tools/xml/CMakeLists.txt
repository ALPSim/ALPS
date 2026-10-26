if(UNIX)
  file(RELATIVE_PATH ALPS_XML_FROM_BINDIR
    "${CMAKE_INSTALL_FULL_BINDIR}" "${CMAKE_INSTALL_FULL_DATADIR}/alps/xml")
  configure_file(alps-xml.in alps-xml @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/alps-xml"
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)

  # Historical names remain aliases of alps-xml.
  foreach(alias IN ITEMS
      "plot2text:plot text" "plot2html:plot html" "plot2gp:plot gnuplot"
      "plot2mpl:plot matplotlib" "plot2xmgr:plot grace"
      "convert2text:convert text" "convert2html:convert html"
      "extracttext:extract text" "extracthtml:extract html"
      "extractgp:extract gnuplot" "extractmpl:extract matplotlib"
      "extractxmgr:extract grace")
    string(REPLACE ":" ";" alias "${alias}")
    list(GET alias 0 ALPS_XML_ALIAS)
    list(GET alias 1 ALPS_XML_COMMAND)
    configure_file(legacy-alias.in ${ALPS_XML_ALIAS} @ONLY NEWLINE_STYLE UNIX)
    install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${ALPS_XML_ALIAS}"
      DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
  endforeach()
  # Deprecated helpers without an alps-xml equivalent warn but keep working.
  configure_file(use_local_stylesheet use_local_stylesheet @ONLY NEWLINE_STYLE UNIX)
  install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/use_local_stylesheet"
    archivecat xslttransform
    DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT xml)
endif()
