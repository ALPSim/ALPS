#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if(LAPACK_FOUND)
  add_definitions(${LAPACK_DEFINITIONS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  add_executable(dwa dwa.cpp)
  target_link_libraries(dwa alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
  install(TARGETS dwa RUNTIME DESTINATION bin COMPONENT applications)
else(LAPACK_FOUND)
  message(STATUS "dwa will not be built since lapack library has not been found")
endif(LAPACK_FOUND)

      IF(ALPS_BUILD_PYTHON)
        add_library (worldlines_c.so SHARED worldlines.cpp)
        add_library (bandstructure_c.so SHARED bandstructure.cpp)
        SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
        set_target_properties(worldlines_c.so PROPERTIES PREFIX "")
        set_target_properties(worldlines_c.so PROPERTIES SUFFIX "")
        set_target_properties(bandstructure_c.so PROPERTIES PREFIX "")
        set_target_properties(bandstructure_c.so PROPERTIES SUFFIX "")
        target_link_libraries(worldlines_c.so alps ${ALPS_LIBRARIES})
        target_link_libraries(bandstructure_c.so alps ${ALPS_LIBRARIES}) 
      ENDIF(ALPS_BUILD_PYTHON)
    install(TARGETS worldlines_c.so bandstructure_c.so COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps/dwa
            LIBRARY DESTINATION lib/pyalps/dwa)

set(DWA_LIB_TOOLS
scripts/tools.py
)

     install(FILES ${DWA_LIB_TOOLS} COMPONENT applications DESTINATION lib/pyalps/dwa)

if(FALSE)
if (ALPS_BUILD_PYTHON)
  python_add_module(worldlines worldlines.cpp)
  if(BUILD_SHARED_LIBS)
    set_target_properties(worldlines PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
  if(WIN32 AND NOT UNIX)
    set_target_properties(worldlines PROPERTIES SUFFIX ".pyd")
  endif(WIN32 AND NOT UNIX)
  endif (BUILD_SHARED_LIBS)

  set_target_properties(worldlines PROPERTIES PREFIX "")
  target_link_libraries(worldlines ${ALPS_LIBRARIES} ${BOOST_LINK_LIBS})

  if(ALPS_INSTALL_VISTRAILS_PACKAGES)
      # TODO: check with Matthias
      install(TARGETS worldlines COMPONENT vistrails
              RUNTIME DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}
              ARCHIVE DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}
              LIBRARY DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR})
  endif(ALPS_INSTALL_VISTRAILS_PACKAGES)

  if(NOT ALPS_FOR_VISTRAILS)
    install(TARGETS worldlines COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib/pyalps/dwa
            LIBRARY DESTINATION lib/pyalps/dwa)
  endif(NOT ALPS_FOR_VISTRAILS)
endif(ALPS_BUILD_PYTHON)
endif(FALSE)
