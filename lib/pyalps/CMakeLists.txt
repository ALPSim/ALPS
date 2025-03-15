#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

#
# python exports
#

set(ALPS_SHARED_CPPFLAGS PYALPS_EXPORTS=1)
set(ALPS_STATIC_CPPFLAGS "")

if (ALPS_HAVE_PYTHON AND NOT ALPS_BUILD_LIBS_ONLY)

  set(OLD_SHARED ${BUILD_SHARED_LIBS})
  set(BUILD_SHARED_LIBS ON)
  set(PYALEA_SOURCES ../../src/alps/python/pyalea.cpp )
  set(PYMCDATA_SOURCES ../../src/alps/python/pymcdata.cpp )
  set(PYTOOLS_SOURCES ../../src/alps/python/pytools.cpp)

  set(PYALPS_SOURCES pyalea_c pymcdata_c pytools_c pyngsparams_c pyngshdf5_c pyngsbase_c  
      pyngsobservable_c pyngsobservables_c pyngsresult_c pyngsresults_c pyngsapi_c pyngsrandom01_c
  )

  if(ALPS_NGS_USE_NEW_ALEA)
    list(APPEND PYALPS_SOURCES pyngsaccumulator_c)
  endif(ALPS_NGS_USE_NEW_ALEA)

  set (MAXENT_SOURCES ../../tool/maxent.cpp ../../tool/maxent_helper.cpp
                      ../../tool/maxent_simulation.cpp ../../tool/maxent_parms.cpp)

  set (CTHYB_SOURCES ../../applications/dmft/qmc/hybridization/hybmain.cpp
       ../../applications/dmft/qmc/hybridization/hybsim.cpp
       ../../applications/dmft/qmc/hybridization/hyblocal.cpp
       ../../applications/dmft/qmc/hybridization/hybint.cpp
       ../../applications/dmft/qmc/hybridization/hybfun.cpp
       ../../applications/dmft/qmc/hybridization/hybretintfun.cpp
       ../../applications/dmft/qmc/hybridization/hybmatrix.cpp
       ../../applications/dmft/qmc/hybridization/hybmatrix_ft.cpp
       ../../applications/dmft/qmc/hybridization/hybconfig.cpp
       ../../applications/dmft/qmc/hybridization/hybupdates.cpp
       ../../applications/dmft/qmc/hybridization/hybevaluate.cpp
       ../../applications/dmft/qmc/hybridization/hybmeasurements.cpp)
  set (CTINT_SOURCES ../../applications/dmft/qmc/interaction_expansion2/main.cpp
       ../../applications/dmft/qmc/fouriertransform.C
       ../../applications/dmft/qmc/interaction_expansion2/auxiliary.cpp
       ../../applications/dmft/qmc/interaction_expansion2/observables.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/fastupdate.cpp
       ../../applications/dmft/qmc/interaction_expansion2/selfenergy.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/solver.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/io.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/splines.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/interaction_expansion.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/measurements.cpp 
       ../../applications/dmft/qmc/interaction_expansion2/model.cpp)

  set(PYNGSPARAMS_SOURCES ../../src/alps/ngs/python/params.cpp)
  set(PYNGSHDF5_SOURCES ../../src/alps/ngs/python/hdf5.cpp)
  set(PYNGSBASE_SOURCES ../../src/alps/ngs/python/mcbase.cpp)
  set(PYNGSOBSERVABLE_SOURCES ../../src/alps/ngs/python/observable.cpp)
  set(PYNGSOBSERVABLES_SOURCES ../../src/alps/ngs/python/observables.cpp)
  set(PYNGSRESULT_SOURCES ../../src/alps/ngs/python/result.cpp)
  set(PYNGSRESULTS_SOURCES ../../src/alps/ngs/python/results.cpp)
  set(PYNGSAPI_SOURCES ../../src/alps/ngs/python/api.cpp)
  set(PYNGSRANDOM01_SOURCES ../../src/alps/ngs/python/random01.cpp)
  set(PYNGSACCUMULATOR_SOURCES ../../src/alps/ngs/python/accumulator.cpp)

  if(LAPACK_FOUND AND ALPS_BUILD_APPLICATIONS)
    set(PYALPS_SOURCES ${PYALPS_SOURCES} maxent_c dwa_c cthyb ctint)
    python_add_module(maxent_c ${MAXENT_SOURCES})
    python_add_module(cthyb ${CTHYB_SOURCES})
    python_add_module(ctint ${CTINT_SOURCES})
    python_add_module(dwa_c ../../applications/qmc/dwa/python/dwa.cpp)
    include_directories(../../applications/qmc/dwa)
    include_directories(../../applications/dmft/qmc)
    set_target_properties(maxent_c PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE") 
    set_target_properties(cthyb PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE") 
    set_target_properties(ctint PROPERTIES COMPILE_FLAGS "-DBUILD_PYTHON_MODULE") 
  endif(LAPACK_FOUND AND ALPS_BUILD_APPLICATIONS)

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

  if(ALPS_NGS_USE_NEW_ALEA)
    python_add_module(pyngsaccumulator_c ${PYNGSACCUMULATOR_SOURCES})
  endif(ALPS_NGS_USE_NEW_ALEA)



  FOREACH (name ${PYALPS_SOURCES})
    if(BUILD_SHARED_LIBS)
      set_target_properties(${name} PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
      set_target_properties(${name} PROPERTIES PREFIX "")
      if(WIN32 AND NOT UNIX)
        set_target_properties(${name} PROPERTIES SUFFIX ".pyd")
      endif(WIN32 AND NOT UNIX)
    endif(BUILD_SHARED_LIBS)
    target_link_libraries(${name} ${LINK_LIBRARIES} ${BLAS_LIBRARY} ${LAPACK_LIBRARY} ${LAPACK_LINKER_FLAGS})
    if(ALPS_PYTHON_WHEEL)
      target_link_libraries(${name} alps_python)
      set_target_properties(${name} PROPERTIES INSTALL_RPATH "$ORIGIN/lib" )
    else()
      target_link_libraries(${name} alps)
      set_target_properties(${name} PROPERTIES INSTALL_RPATH "$ORIGIN/../../.." )
    endif()
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
    if(NOT ALPS_PYTHON_WHEEL)
      install(TARGETS ${PYALPS_SOURCES}
            COMPONENT python
            RUNTIME DESTINATION ${ALPS_PYTHON_LIB_DEST_ROOT}/pyalps/bin
            ARCHIVE DESTINATION ${ALPS_PYTHON_LIB_DEST_ROOT}/pyalps
            LIBRARY DESTINATION ${ALPS_PYTHON_LIB_DEST_ROOT}/pyalps)
    else(NOT ALPS_PYTHON_WHEEL)
      install(TARGETS ${PYALPS_SOURCES}
            COMPONENT python
            RUNTIME DESTINATION pyalps/bin
            ARCHIVE DESTINATION pyalps
            LIBRARY DESTINATION pyalps)
    endif(NOT ALPS_PYTHON_WHEEL)
  endif(NOT ALPS_FOR_VISTRAILS)
  set(BUILD_SHARED_LIBS ${OLD_SHARED})
endif (ALPS_HAVE_PYTHON AND NOT ALPS_BUILD_LIBS_ONLY)
