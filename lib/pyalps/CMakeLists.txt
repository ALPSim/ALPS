#
# pyalea
#

if (ALPS_HAVE_PYTHON)
set(PYALEA_SOURCES "value_with_error.C")
endif (ALPS_HAVE_PYTHON)


if(BUILD_SHARED_LIBS)
  add_library(pyalea ${PYALEA_SOURCES})
  set_target_properties(pyalea PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
  set_target_properties(pyalea PROPERTIES  PREFIX "" SUFFIX ".so")
endif(BUILD_SHARED_LIBS)

if(ALPS_FOR_VISTRAILS)
  target_link_libraries(pyalea boost ${LINK_LIBRARIES})
else(ALPS_FOR_VISTRAILS)
  target_link_libraries(pyalea alps boost ${LINK_LIBRARIES})
endif(ALPS_FOR_VISTRAILS)

#######################################################################
# install
#######################################################################
install(TARGETS pyalea
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib/pyalps
        LIBRARY DESTINATION lib/pyalps)

