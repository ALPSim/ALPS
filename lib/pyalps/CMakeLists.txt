#
# pyalea
#

if (ALPS_HAVE_PYTHON)
  set(PYALEA_SOURCES value_with_error.C)
  set(PYTOOLS_SOURCES pytools.C)


  if(BUILD_SHARED_LIBS)
    python_add_module(pyalea ${PYALEA_SOURCES})
    python_add_module(pytools ${PYTOOLS_SOURCES})
    set_target_properties(pyalea PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pytools PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pyalea PROPERTIES PREFIX "")
    set_target_properties(pytools PROPERTIES PREFIX "")
  if(WIN32 AND NOT UNIX)
    set_target_properties(pyalea PROPERTIES SUFFIX ".pyd")
    set_target_properties(pytools PROPERTIES SUFFIX ".pyd")
  endif(WIN32 AND NOT UNIX)
  endif(BUILD_SHARED_LIBS)

  target_link_libraries(pyalea boost ${LINK_LIBRARIES})
  target_link_libraries(pytools alps boost ${LINK_LIBRARIES})

  #######################################################################
  # install
  #######################################################################
  if(ALPS_FOR_VISTRAILS)
    install(TARGETS pyalea pytools
            RUNTIME DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            ARCHIVE DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  else(ALPS_FOR_VISTRAILS)
    install(TARGETS pyalea pytools
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps
            LIBRARY DESTINATION lib/pyalps)
  endif(ALPS_FOR_VISTRAILS)

endif (ALPS_HAVE_PYTHON)
