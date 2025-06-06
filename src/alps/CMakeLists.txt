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
# libalps
#

set(ALPS_SOURCES "")
set(ALPS_PYTHON_SOURCES "")
set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
set(ALPS_STATIC_CPPFLAGS "")

set(ALPS_SOURCES
    # utility
    utility/encode.cpp utility/copyright.cpp utility/os.cpp utility/temporary_filename.cpp
    utility/vmusage.cpp
    # alea
    alea/nan.C alea/observable.C alea/observablefactory.C alea/observableset.C
    # expression
    expression/evaluator.C
    # lattice
    lattice/disorder.C lattice/latticedescriptor.C lattice/latticegraphdescriptor.C
    lattice/latticelibrary.C lattice/unitcell.C
    # model
    model/bondterm.C model/globaloperator.C model/modellibrary.C model/siteterm.C
    # osiris
    osiris/buffer.C osiris/dump.C osiris/dumparchive.C osiris/process.C osiris/comm.C
    osiris/mpdump.C osiris/xdrdump.C
    # parameter
    parameter/parameter.C parameter/parameterlist.C parameter/parameters.C
    # parser
    parser/parser.C parser/xmlattributes.C parser/xmlhandler.C parser/xmlparser.C
    parser/xmlstream.C parser/xslt_path.C
    #hdf5
    hdf5/archive.cpp ngs/lib/signal.cpp ngs/lib/stacktrace.cpp
    # random
    random/rngfactory.C random/parallel/get_prime.cpp random/parallel/primelist_64.cpp
    # scheduler
    scheduler/abstract_task.C scheduler/convertxdr.C scheduler/factory.C scheduler/info.C
    scheduler/master_scheduler.C scheduler/measurement_operators.C
    scheduler/montecarlo.C scheduler/mpp_scheduler.C scheduler/options.C
    scheduler/remote_task.C scheduler/remote_worker.C scheduler/scheduler.C
    scheduler/serial_scheduler.C scheduler/signal.C scheduler/single_scheduler.C
    scheduler/slave_task.C scheduler/task.C scheduler/worker.C scheduler/workertask.C
    # ngs sources
    stop_callback.cpp ngs/lib/api.cpp ngs/lib/short_print.cpp ngs/lib/sleep.cpp
    # alea
    # parameters
    ngs/lib/paramvalue.cpp ngs/lib/paramproxy.cpp ngs/lib/params.cpp
    # scheduler
    mcbase.cpp parseargs.cpp
#deprecated, TODO: replace by new alea
    ngs/lib/mcoptions.cpp ngs/lib/observablewrappers.cpp ngs/lib/make_deprecated_parameters.cpp ngs/lib/make_parameters_from_xml.cpp
    ngs/lib/mcresult.cpp ngs/lib/mcresults.cpp ngs/lib/mcobservables.cpp ngs/lib/mcobservable.cpp
  )

  if (NOT ALPS_HAVE_RPC_XDR_H)
    set(ALPS_SOURCES ${ALPS_SOURCES} osiris/xdr.c osiris/xdr_array.c osiris/xdr_float.c osiris/xdr_stdio.c)
  endif (NOT ALPS_HAVE_RPC_XDR_H)
  
  if (ALPS_BUILD_PYTHON)
    set(ALPS_PYTHON_SOURCES ${ALPS_PYTHON_SOURCES} ngs/lib/get_numpy_type.cpp hdf5/python.cpp python/numpy_array.cpp)
  endif (ALPS_BUILD_PYTHON)

  # OpenMPI ULFM
  if (ALPS_NGS_OPENMPI_ULFM)
    set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/ulfm.cpp)
  endif (ALPS_NGS_OPENMPI_ULFM)

  # new alea
  if (ALPS_NGS_USE_NEW_ALEA)
#    set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/alea/accumulator_set.cpp ngs/lib/alea/result_set.cpp)
    set(ALPS_SOURCES ${ALPS_SOURCES} ngs/accumulator/accumulator.cpp)
  endif (ALPS_NGS_USE_NEW_ALEA)

if (NOT ALPS_FOR_VISTRAILS)
  set(ALPS_SOURCES ${ALPS_SOURCES}
      # parapack
      parapack/clone.C parapack/clone_info.C parapack/filelock.C parapack/logger.C
      parapack/job.C parapack/mc_worker.C parapack/measurement.C parapack/process_impl.C
      parapack/option.C parapack/parapack.C parapack/queue.C parapack/rng_helper.C
      parapack/types.C parapack/util.C parapack/version.C parapack/worker_factory.C
      parapack/staging.C
      # ngs-parapack
      ngs/lib/clone.cpp ngs/lib/clone_info.cpp ngs/lib/job.cpp ngs/lib/parapack.cpp ngs/lib/worker_factory.cpp
  )
endif(NOT ALPS_FOR_VISTRAILS)

if(ALPS_HAVE_PYTHON)
  if(ALPS_PYTHON_WHEEL)
    set(ALPS_PYTHON_SOURCES ${ALPS_PYTHON_SOURCES} ${ALPS_SOURCES})
    add_library(alps_python ${ALPS_PYTHON_SOURCES})
    target_compile_definitions(alps_python PRIVATE ALPS_HAVE_PYTHON)
    get_target_property(XXX alps_python COMPILE_DEFINITIONS)
    message(STATUS "ALPS_PYTHON: ${XXX}")
    if(ALPS_HAVE_BOOST_NUMPY)
      target_compile_definitions(alps_python INTERFACE ALPS_HAVE_BOOST_NUMPY)
    endif()
  else()
    set(ALPS_SOURCES ${ALPS_PYTHON_SOURCES} ${ALPS_SOURCES})
  endif()
endif()

add_library(alps ${ALPS_SOURCES})

find_package(HDF5)
if (Boost_FOUND)
  # link to ${Boost_LIBRARIES} when precompiled Boost libraries found
  set(ALPS_LINK_LIBS ${Boost_LIBRARIES} ${HDF5_LIBRARIES}) # ${SZIP_LIBRARIES} ${ZLIB_LIBRARIES})
  if(PYTHONLIBS_FOUND)
    if(ALPS_PYTHON_WHEEL AND PYTHON_VERSION VERSION_GREATER_EQUAL "3.13" AND APPLE)
      set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${PYTHON_LIBRARY}) # ${PYTHON_EXTRA_LIBS})
    else(ALPS_PYTHON_WHEEL AND PYTHON_VERSION VERSION_GREATER_EQUAL "3.13" AND APPLE)
      set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${PYTHON_LIBRARY} ${PYTHON_EXTRA_LIBS})
    endif(ALPS_PYTHON_WHEEL AND PYTHON_VERSION VERSION_GREATER_EQUAL "3.13" AND APPLE)
  endif(PYTHONLIBS_FOUND)
  if(MPI_FOUND)
    set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${MPI_LIBRARIES})
    if(MPI_EXTRA_LIBRARY)
      set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${MPI_EXTRA_LIBRARY})
    endif(MPI_EXTRA_LIBRARY)
  endif(MPI_FOUND)
  target_link_libraries(alps ${ALPS_LINK_LIBS})
else (Boost_FOUND)
  # "boost" target available when Boost libraries are built from source
  target_link_libraries(alps ${ALPS_BOOST_LIBRARY_NAME} ${HDF5_LIBRARIES}) # ${SZIP_LIBRARIES} ${ZLIB_LIBRARIES})
  if(ALPS_HAVE_PYTHON)
    if(ALPS_PYTHON_WHEEL)
      target_link_libraries(alps_python ${ALPS_BOOST_LIBRARY_NAME} ${ALPS_BOOST_PYTHON_LIBRARY_NAME} ${HDF5_LIBRARIES}) # ${SZIP_LIBRARIES} ${ZLIB_LIBRARIES})
    else()
      target_link_libraries(alps ${ALPS_BOOST_PYTHON_LIBRARY_NAME}) # ${SZIP_LIBRARIES} ${ZLIB_LIBRARIES})
    endif()
  endif()
endif (Boost_FOUND)

if(BUILD_SHARED_LIBS)
  if(HDF5_DEFINITIONS)
    foreach (_hdf5_def ${HDF5_DEFINITIONS})
      string(REPLACE "-D" "" _escaped_def ${_hdf5_def})
      list(APPEND ALPS_SHARED_CPPFLAGS ${_escaped_def})
    endforeach(_hdf5_def)
  endif(HDF5_DEFINITIONS)
  
  set_target_properties(alps PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
  if(ALPS_HAVE_PYTHON)
    if(ALPS_PYTHON_WHEEL)
      set_target_properties(alps_python PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
      target_compile_definitions(alps_python PUBLIC ALPS_HAVE_PYTHON)
      if(ALPS_HAVE_BOOST_NUMPY)
        target_compile_definitions(alps_python PUBLIC ALPS_HAVE_BOOST_NUMPY)
      endif()
    else(ALPS_PYTHON_WHEEL)
      target_compile_definitions(alps PUBLIC ALPS_HAVE_PYTHON)
      if(ALPS_HAVE_BOOST_NUMPY)
        target_compile_definitions(alps PUBLIC ALPS_HAVE_BOOST_NUMPY)
      endif()
    endif(ALPS_PYTHON_WHEEL)
  endif()
else(BUILD_SHARED_LIBS)
  set_target_properties(alps PROPERTIES COMPILE_DEFINITIONS "${ALPS_STATIC_CPPFLAGS}")
  if(ALPS_HAVE_PYTHON)
    set_target_properties(alps_python PROPERTIES COMPILE_DEFINITIONS "${ALPS_STATIC_CPPFLAGS}")
  endif()
endif(BUILD_SHARED_LIBS)


if(MSVC)
  get_target_property(LIB_NAME alps LOCATION)
  add_custom_command(TARGET alps POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${LIB_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)

# Set soversion for library
if(NOT WIN32 AND NOT APPLE AND NOT ALPS_PYTHON_WHEEL)
  set_target_properties(alps
    PROPERTIES
    SOVERSION "${ALPS_VERSION_MAJOR}"
    VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}" )
endif(NOT WIN32 AND NOT APPLE AND NOT ALPS_PYTHON_WHEEL)


#boost librt linking
if(NOT WIN32 AND NOT APPLE)
    set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} rt)
endif(NOT WIN32 AND NOT APPLE)

if(ALPS_BUILD_FORTRAN)
  set(ALPS_FORTRAN_SOURCES fortran/fwrapper_impl.C)
  add_library(alps_fortran STATIC ${ALPS_FORTRAN_SOURCES})
  target_link_libraries(alps_fortran ${CMAKE_Fortran_RUNTIME_LIBRARIES})
  set_target_properties(alps_fortran PROPERTIES COMPILE_DEFINITIONS "${ALPS_STATIC_CPPFLAGS}")
endif()
#######################################################################
# install
#######################################################################
if(NOT ALPS_PYTHON_WHEEL)
if (NOT ALPS_FOR_VISTRAILS)
  install(TARGETS alps COMPONENT libraries
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION bin)

endif (NOT ALPS_FOR_VISTRAILS)

if (ALPS_FOR_VISTRAILS)
  if (APPLE)
   SET_TARGET_PROPERTIES(alps
    PROPERTIES BUILD_WITH_INSTALL_RPATH ON
     INSTALL_NAME_DIR "@executable_path/../Frameworks"
    )
  endif(APPLE)
  message (STATUS "Putting into ${VISTRAILS_DYLIB_DIR}")
  install(TARGETS alps  COMPONENT vistrails
          ARCHIVE DESTINATION ${VISTRAILS_DYLIB_DIR}
          LIBRARY DESTINATION ${VISTRAILS_DYLIB_DIR}
          RUNTIME DESTINATION ${VISTRAILS_DYLIB_DIR})
endif (ALPS_FOR_VISTRAILS)


if(ALPS_BUILD_FORTRAN)
  if(NOT ALPS_FOR_VISTRAILS)
    install(TARGETS alps_fortran COMPONENT libraries
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
            RUNTIME DESTINATION bin)
  endif(NOT ALPS_FOR_VISTRAILS)
  if(ALPS_INSTALL_VISTRAILS_PACKAGES)
    install(TARGETS alps_fortran COMPONENT vistrails
            ARCHIVE DESTINATION ${VISTRAILS_DYLIB_DIR}
            LIBRARY DESTINATION ${VISTRAILS_DYLIB_DIR}
            RUNTIME DESTINATION ${VISTRAILS_DYLIB_DIR})
  endif(ALPS_INSTALL_VISTRAILS_PACKAGES)
endif(ALPS_BUILD_FORTRAN)
else ()
  install(TARGETS alps_python COMPONENT libraries
          ARCHIVE DESTINATION pyalps/lib
          LIBRARY DESTINATION pyalps/lib
          RUNTIME DESTINATION pyalps/bin)
  install(TARGETS alps COMPONENT libraries
          ARCHIVE DESTINATION pyalps/lib
          LIBRARY DESTINATION pyalps/lib
          RUNTIME DESTINATION pyalps/bin)
endif()
