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


  python_add_module(pyalea_c ${PYALEA_SOURCES})
  python_add_module(pytools_c ${PYTOOLS_SOURCES})
  if(BUILD_SHARED_LIBS)
    set_target_properties(pyalea_c PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pytools_c PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
    set_target_properties(pyalea_c PROPERTIES PREFIX "")
    set_target_properties(pytools_c PROPERTIES PREFIX "")
    if(WIN32 AND NOT UNIX)
      set_target_properties(pyalea_c PROPERTIES SUFFIX ".pyd")
      set_target_properties(pytools_c PROPERTIES SUFFIX ".pyd")
    endif(WIN32 AND NOT UNIX)
  endif(BUILD_SHARED_LIBS)
  target_link_libraries(pyalea_c alps ${LINK_LIBRARIES})
  target_link_libraries(pytools_c alps ${LINK_LIBRARIES})


  #######################################################################
  # install
  #######################################################################
  if (ALPS_INSTALL_VISTRAILS_PACKAGES)
    install(TARGETS pyalea_c pytools_c COMPONENT vistrails
            RUNTIME DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            ARCHIVE DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  endif (ALPS_INSTALL_VISTRAILS_PACKAGES)
  if (NOT ALPS_FOR_VISTRAILS)
    install(TARGETS pyalea_c pytools_c COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps
            LIBRARY DESTINATION lib/pyalps)
  endif(NOT ALPS_FOR_VISTRAILS)

endif (ALPS_HAVE_PYTHON)
