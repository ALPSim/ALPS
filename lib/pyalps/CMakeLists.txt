#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# pyalea
#

if (ALPS_HAVE_PYTHON)
  set(PYALEA_SOURCES ../../src/alps/python/pyalea.cpp ../../src/alps/python/numpy_array.cpp)
  set(PYMCDATA_SOURCES ../../src/alps/python/pymcdata.cpp ../../src/alps/python/numpy_array.cpp)
  set(PYTOOLS_SOURCES ../../src/alps/python/pytools.cpp)
  set(PYHDF5_SOURCES ../../src/alps/python/pyhdf5.cpp ../../src/alps/python/numpy_array.cpp ../../src/alps/hdf5/hdf5.cpp)

  python_add_module(pyalea_c ${PYALEA_SOURCES})
  python_add_module(pymcdata_c ${PYMCDATA_SOURCES})
  python_add_module(pytools_c ${PYTOOLS_SOURCES})
  python_add_module(pyhdf5_c ${PYHDF5_SOURCES})

  FOREACH (name pyhdf5_c)
    if(BUILD_SHARED_LIBS)
      set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
      set_target_properties(${name} PROPERTIES PREFIX "")
      if(WIN32 AND NOT UNIX)
        set_target_properties(${name} PROPERTIES SUFFIX ".pyd")
      endif(WIN32 AND NOT UNIX)
    endif(BUILD_SHARED_LIBS)
    find_package(HDF5)
    target_link_libraries(${name} boost ${LINK_LIBRARIES} ${HDF5_LIBRARIES})
  ENDFOREACH(name)


  FOREACH (name pyalea_c pymcdata_c pytools_c)
    if(BUILD_SHARED_LIBS)
      set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
      set_target_properties(${name} PROPERTIES PREFIX "")
      if(WIN32 AND NOT UNIX)
        set_target_properties(${name} PROPERTIES SUFFIX ".pyd")
      endif(WIN32 AND NOT UNIX)
    endif(BUILD_SHARED_LIBS)
    target_link_libraries(${name} alps ${LINK_LIBRARIES})
  ENDFOREACH(name)

  #######################################################################
  # install
  #######################################################################
  if (ALPS_INSTALL_VISTRAILS_PACKAGES)
    install(TARGETS pyalea_c pymcdata_c pytools_c pyhdf5_c COMPONENT vistrails
            RUNTIME DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            ARCHIVE DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  endif (ALPS_INSTALL_VISTRAILS_PACKAGES)
  if (NOT ALPS_FOR_VISTRAILS)
    install(TARGETS pyalea_c pymcdata_c pytools_c pyhdf5_c COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps
            LIBRARY DESTINATION lib/pyalps)
  endif(NOT ALPS_FOR_VISTRAILS)

endif (ALPS_HAVE_PYTHON)
