#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# python exports
#

set(ALPS_SHARED_CPPFLAGS PYALPS_EXPORTS=1)
set(ALPS_STATIC_CPPFLAGS "")

if (ALPS_HAVE_PYTHON)

  set(PYALEA_SOURCES ../../src/alps/python/pyalea.cpp )
  set(PYMCDATA_SOURCES ../../src/alps/python/pymcdata.cpp )
  set(PYTOOLS_SOURCES ../../src/alps/python/pytools.cpp)

  set(PYALPS_SOURCES pyalea_c pymcdata_c pytools_c pyngsparams_c pyngshdf5_c pyngsbase_c  
      pyngsobservable_c pyngsobservables_c pyngsresult_c pyngsresults_c pyngsapi_c pyngsrandom01_c
  )

  set (MAXENT_SOURCES ../../tool/maxent.cpp ../../tool/maxent_helper.cpp
                      ../../tool/maxent_simulation.cpp ../../tool/maxent_parms.cpp)




  set(PYNGSPARAMS_SOURCES ../../src/alps/ngs/python/params.cpp)
  set(PYNGSHDF5_SOURCES ../../src/alps/ngs/python/hdf5.cpp)
  set(PYNGSBASE_SOURCES ../../src/alps/ngs/python/mcbase.cpp)
  set(PYNGSOBSERVABLE_SOURCES ../../src/alps/ngs/python/observable.cpp)
  set(PYNGSOBSERVABLES_SOURCES ../../src/alps/ngs/python/observables.cpp)
  set(PYNGSRESULT_SOURCES ../../src/alps/ngs/python/result.cpp)
  set(PYNGSRESULTS_SOURCES ../../src/alps/ngs/python/results.cpp)
  set(PYNGSAPI_SOURCES ../../src/alps/ngs/python/api.cpp)
  set(PYNGSRANDOM01_SOURCES ../../src/alps/ngs/python/random01.cpp)


  if(LAPACK_FOUND)
    set(PYALPS_SOURCES ${PYALPS_SOURCES} maxent_c dwa_c)
    python_add_module(maxent_c ${MAXENT_SOURCES})
    python_add_module(dwa_c ../../applications/qmc/dwa/python/dwa.cpp)
    include_directories(../../applications/qmc/dwa)
    set_target_properties(maxent_c PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE") 
 endif(LAPACK_FOUND)

  python_add_module(pyalea_c ${PYALEA_SOURCES})
  python_add_module(pymcdata_c ${PYMCDATA_SOURCES})
  python_add_module(pytools_c ${PYTOOLS_SOURCES})
  python_add_module(pyngsparams_c ${PYNGSPARAMS_SOURCES})
  python_add_module(pyngshdf5_c ${PYNGSHDF5_SOURCES})
  python_add_module(pyngsbase_c ${PYNGSBASE_SOURCES})
  python_add_module(pyngsobservable_c ${PYNGSOBSERVABLE_SOURCES})
  python_add_module(pyngsobservables_c ${PYNGSOBSERVABLES_SOURCES})
  python_add_module(pyngsresult_c ${PYNGSRESULT_SOURCES})
  python_add_module(pyngsresults_c ${PYNGSRESULTS_SOURCES})
  python_add_module(pyngsapi_c ${PYNGSAPI_SOURCES})
  python_add_module(pyngsrandom01_c ${PYNGSRANDOM01_SOURCES})

  FOREACH (name ${PYALPS_SOURCES})
    if(BUILD_SHARED_LIBS)
      set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
      set_target_properties(${name} PROPERTIES PREFIX "")
      if(WIN32 AND NOT UNIX)
        set_target_properties(${name} PROPERTIES SUFFIX ".pyd")
      endif(WIN32 AND NOT UNIX)
    endif(BUILD_SHARED_LIBS)
    target_link_libraries(${name} alps ${LINK_LIBRARIES} ${BLAS_LIBRARY} ${LAPACK_LIBRARY} ${LAPACK_LINKER_FLAGS})
  ENDFOREACH(name)


  #######################################################################
  # install
  #######################################################################
  if (ALPS_INSTALL_VISTRAILS_PACKAGES)
    install(TARGETS ${PYALPS_SOURCES}
            COMPONENT vistrails
            RUNTIME DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            ARCHIVE DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  endif (ALPS_INSTALL_VISTRAILS_PACKAGES)
  if (NOT ALPS_FOR_VISTRAILS)
    install(TARGETS ${PYALPS_SOURCES}
            COMPONENT python
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION ${PYTHON_EXTENSION_DIR}/pyalps
            LIBRARY DESTINATION ${PYTHON_EXTENSION_DIR}/pyalps)
  endif(NOT ALPS_FOR_VISTRAILS)

endif (ALPS_HAVE_PYTHON)
