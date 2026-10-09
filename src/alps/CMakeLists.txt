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
set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
set(ALPS_STATIC_CPPFLAGS "")

set(ALPS_SOURCES
    # utility
    utilities/src/encode.cpp utilities/src/copyright.cpp utilities/src/os.cpp utilities/src/temporary_filename.cpp
    utilities/src/vmusage.cpp
    utilities/src/citations.cpp
    cli/src/cli.cpp
    # alea
    alea/src/alea/nan.C alea/src/alea/observable.C alea/src/alea/observablefactory.C alea/src/alea/observableset.C
    # expression
    expression/src/expression/evaluator.C
    # lattice
    lattice/src/lattice/disorder.C lattice/src/lattice/latticedescriptor.C lattice/src/lattice/latticegraphdescriptor.C
    lattice/src/lattice/latticelibrary.C lattice/src/lattice/unitcell.C
    # model
    model/src/model/bondterm.C model/src/model/globaloperator.C model/src/model/modellibrary.C model/src/model/siteterm.C
    # osiris
    osiris/src/osiris/buffer.C osiris/src/osiris/dump.C osiris/src/osiris/dumparchive.C osiris/src/osiris/process.C osiris/src/osiris/comm.C
    osiris/src/osiris/mpdump.C osiris/src/osiris/xdrdump.C
    # parameter
    legacy_parameters/src/parameter/parameter.C legacy_parameters/src/parameter/parameterlist.C legacy_parameters/src/parameter/parameters.C
    # parser
    xml/src/parser/parser.C xml/src/parser/xmlattributes.C xml/src/parser/xmlhandler.C xml/src/parser/xmlparser.C
    xml/src/parser/xmlstream.C xml/src/parser/xslt_path.C
    #hdf5
    hdf5/src/archive.cpp hdf5/src/signal.cpp utilities/src/stacktrace.cpp
    # random
    random/src/random/rngfactory.C random/src/random/parallel/get_prime.cpp random/src/random/parallel/primelist_64.cpp
    # scheduler
    scheduler/src/scheduler/abstract_task.C scheduler/src/scheduler/convertxdr.C scheduler/src/scheduler/factory.C scheduler/src/scheduler/info.C
    scheduler/src/scheduler/master_scheduler.C scheduler/src/scheduler/measurement_operators.C
    scheduler/src/scheduler/montecarlo.C scheduler/src/scheduler/mpp_scheduler.C scheduler/src/scheduler/options.C
    scheduler/src/scheduler/remote_task.C scheduler/src/scheduler/remote_worker.C scheduler/src/scheduler/scheduler.C
    scheduler/src/scheduler/serial_scheduler.C scheduler/src/scheduler/signal.C scheduler/src/scheduler/single_scheduler.C
    scheduler/src/scheduler/slave_task.C scheduler/src/scheduler/task.C scheduler/src/scheduler/worker.C scheduler/src/scheduler/workertask.C
    # ngs sources
    mc/src/stop_callback.cpp mc/src/ngs/api.cpp utilities/src/short_print.cpp utilities/src/sleep.cpp
    # alea
    # parameters
    params/src/paramvalue.cpp params/src/paramproxy.cpp params/src/params.cpp
    # scheduler
    mc/src/mcbase.cpp cli/src/parseargs.cpp
#deprecated, TODO: replace by new alea
    cli/src/mcoptions.cpp alea/src/ngs/observablewrappers.cpp params/adapters/make_deprecated_parameters.cpp params/adapters/make_parameters_from_xml.cpp
    alea/src/ngs/mcresult.cpp alea/src/ngs/mcresults.cpp alea/src/ngs/mcobservables.cpp alea/src/ngs/mcobservable.cpp
  )

  if (NOT ALPS_HAVE_RPC_XDR_H)
    set(ALPS_SOURCES ${ALPS_SOURCES} ../../third_party/xdr/src/xdr.c ../../third_party/xdr/src/xdr_array.c ../../third_party/xdr/src/xdr_float.c ../../third_party/xdr/src/xdr_stdio.c)
  endif (NOT ALPS_HAVE_RPC_XDR_H)
  
  # OpenMPI ULFM
  if (ALPS_NGS_OPENMPI_ULFM)
    set(ALPS_SOURCES ${ALPS_SOURCES} scheduler/src/ngs/ulfm.cpp)
  endif (ALPS_NGS_OPENMPI_ULFM)

  # new alea
  if (ALPS_NGS_USE_NEW_ALEA)
#    set(ALPS_SOURCES ${ALPS_SOURCES} ngs/lib/alea/accumulator_set.cpp ngs/lib/alea/result_set.cpp)
    set(ALPS_SOURCES ${ALPS_SOURCES} compatibility/ngs/accumulator/accumulator.cpp)
  endif (ALPS_NGS_USE_NEW_ALEA)

set(ALPS_SOURCES ${ALPS_SOURCES}
    # parapack
    parapack/src/parapack/clone.C parapack/src/parapack/clone_info.C parapack/src/parapack/filelock.C parapack/src/parapack/logger.C
    parapack/src/parapack/job.C parapack/src/parapack/mc_worker.C parapack/src/parapack/measurement.C parapack/src/parapack/process_impl.C
    parapack/src/parapack/option.C parapack/src/parapack/parapack.C parapack/src/parapack/queue.C parapack/src/parapack/rng_helper.C
    parapack/src/parapack/types.C parapack/src/parapack/util.C parapack/src/parapack/version.C parapack/src/parapack/worker_factory.C
    parapack/src/parapack/staging.C
    # ngs-parapack
    parapack/src/ngs/clone.cpp parapack/src/ngs/clone_info.cpp parapack/src/ngs/job.cpp parapack/src/ngs/parapack.cpp parapack/src/ngs/worker_factory.cpp
)

add_library(alps ${ALPS_SOURCES})

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # The legacy "lagged_fibonacci607" generator uses a 48-bit word with a 32-bit
  # result type, so max() truncates by design; keep it for stream compatibility.
  set_source_files_properties(random/src/random/rngfactory.C PROPERTIES
    COMPILE_OPTIONS "-Wno-constant-conversion")
endif()

find_package(HDF5)
if (Boost_FOUND)
  # link to ${Boost_LIBRARIES} when precompiled Boost libraries found
  set(ALPS_LINK_LIBS ${Boost_LIBRARIES} ${HDF5_LIBRARIES}) # ${SZIP_LIBRARIES})
  if(MPI_FOUND)
    set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${MPI_LIBRARIES})
    if(MPI_EXTRA_LIBRARY)
      set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} ${MPI_EXTRA_LIBRARY})
    endif(MPI_EXTRA_LIBRARY)
  endif(MPI_FOUND)
  target_link_libraries(alps ${ALPS_LINK_LIBS})
else (Boost_FOUND)
  # "boost" target available when Boost libraries are built from source
  target_link_libraries(alps ${ALPS_BOOST_LIBRARY_NAME} ${HDF5_LIBRARIES}) # ${SZIP_LIBRARIES})
endif (Boost_FOUND)

if(BUILD_SHARED_LIBS)
  if(HDF5_DEFINITIONS)
    foreach (_hdf5_def ${HDF5_DEFINITIONS})
      string(REPLACE "-D" "" _escaped_def ${_hdf5_def})
      list(APPEND ALPS_SHARED_CPPFLAGS ${_escaped_def})
    endforeach(_hdf5_def)
  endif(HDF5_DEFINITIONS)
  
  set_target_properties(alps PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
else(BUILD_SHARED_LIBS)
  set_target_properties(alps PROPERTIES COMPILE_DEFINITIONS "${ALPS_STATIC_CPPFLAGS}")
endif(BUILD_SHARED_LIBS)


if(MSVC)
  get_target_property(LIB_NAME alps LOCATION)
  add_custom_command(TARGET alps POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${LIB_NAME} ${PROJECT_BINARY_DIR}/bin)
endif(MSVC)

# Set soversion for library
if(NOT WIN32 AND NOT APPLE)
  set_target_properties(alps
    PROPERTIES
    SOVERSION "${ALPS_VERSION_MAJOR}"
    VERSION "${ALPS_VERSION_MAJOR}.${ALPS_VERSION_MINOR}.${ALPS_VERSION_PATCH}" )
endif(NOT WIN32 AND NOT APPLE)


#boost librt linking
if(NOT WIN32 AND NOT APPLE)
    set(ALPS_LINK_LIBS ${ALPS_LINK_LIBS} rt)
endif(NOT WIN32 AND NOT APPLE)

if(ALPS_BUILD_FORTRAN)
  set(ALPS_FORTRAN_SOURCES fortran/src/fwrapper_impl.C)
  add_library(alps_fortran STATIC ${ALPS_FORTRAN_SOURCES})
  target_link_libraries(alps_fortran ${CMAKE_Fortran_RUNTIME_LIBRARIES})
  set_target_properties(alps_fortran PROPERTIES COMPILE_DEFINITIONS "${ALPS_STATIC_CPPFLAGS}")
endif()
#######################################################################
# install
#######################################################################
install(TARGETS alps COMPONENT libraries
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)

if(ALPS_BUILD_FORTRAN)
  install(TARGETS alps_fortran COMPONENT libraries
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib
          RUNTIME DESTINATION bin)
endif(ALPS_BUILD_FORTRAN)
