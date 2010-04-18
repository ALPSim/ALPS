#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# pyalea
#

if (ALPS_HAVE_PYTHON)
  set(PYALEA_SOURCES ../../src/alps/python/pyalea.cpp)
  set(PYTOOLS_SOURCES ../../src/alps/python/pytools.cpp)


  python_add_module(pyalea ${PYALEA_SOURCES})
  python_add_module(pytools ${PYTOOLS_SOURCES})
  if(BUILD_SHARED_LIBS)
    set_target_properties(pyalea PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pytools PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pyalea PROPERTIES PREFIX "")
    set_target_properties(pytools PROPERTIES PREFIX "")
    if(WIN32 AND NOT UNIX)
      set_target_properties(pyalea PROPERTIES SUFFIX ".pyd")
      set_target_properties(pytools PROPERTIES SUFFIX ".pyd")
    endif(WIN32 AND NOT UNIX)
  endif(BUILD_SHARED_LIBS)
  target_link_libraries(pyalea alps ${LINK_LIBRARIES})
  target_link_libraries(pytools alps ${LINK_LIBRARIES})


  #######################################################################
  # install
  #######################################################################
  if (ALPS_INSTALL_VISTRAILS_PACKAGES)
    install(TARGETS pyalea pytools COMPONENT vistrails
            RUNTIME DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            ARCHIVE DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  endif (ALPS_INSTALL_VISTRAILS_PACKAGES)
  if (NOT ALPS_FOR_VISTRAILS)
    install(TARGETS pyalea pytools COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps
            LIBRARY DESTINATION lib/pyalps)
  endif(NOT ALPS_FOR_VISTRAILS)

endif (ALPS_HAVE_PYTHON)
