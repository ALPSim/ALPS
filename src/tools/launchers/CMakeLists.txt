# Preserve the installed environment helpers during the SDK transition.
if(UNIX)
  foreach(name alpsvars.sh alpsvars.csh)
    configure_file(${name}.in ${name} @ONLY)
    install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${name}"
      DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
  endforeach()
else()
  install(PROGRAMS msxsl.exe DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)
endif()
