#
# pyalea
#

if (ALPS_HAVE_PYTHON)
set(PYALEA_SOURCES "value_with_error.C")
endif (ALPS_HAVE_PYTHON)


if(BUILD_SHARED_LIBS)
  python_add_module(pyalea ${PYALEA_SOURCES})
  set_target_properties(pyalea PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
  set_target_properties(pyalea PROPERTIES PREFIX "")
if(WIN32 AND NOT UNIX)
  set_target_properties(pyalea PROPERTIES SUFFIX ".pyd")
endif(WIN32 AND NOT UNIX)
endif(BUILD_SHARED_LIBS)

if(ALPS_FOR_VISTRAILS)
  target_link_libraries(pyalea boost ${LINK_LIBRARIES})
else(ALPS_FOR_VISTRAILS)
  target_link_libraries(pyalea alps boost ${LINK_LIBRARIES})
endif(ALPS_FOR_VISTRAILS)

#######################################################################
# install
#######################################################################
if(ALPS_FOR_VISTRAILS)
  install(TARGETS pyalea
          ARCHIVE DESTINATION ../Resources/lib/python2.5/pyalps
          LIBRARY DESTINATION ../Resources/lib/python2.5/pyalps)
else(ALPS_FOR_VISTRAILS)
  install(TARGETS pyalea
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib/pyalps
          LIBRARY DESTINATION lib/pyalps)
endif(ALPS_FOR_VISTRAILS)

