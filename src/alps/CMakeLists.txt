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

  list(APPEND ALPS_SOURCES osiris/xdr.c osiris/xdr_array.c osiris/xdr_float.c osiris/xdr_stdio.c)
  
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

add_library(alps ${ALPS_SOURCES})

add_library(ALPS::alps ALIAS alps)
add_library(alps_headers INTERFACE)
add_library(ALPS::headers ALIAS alps_headers)
set_target_properties(alps_headers PROPERTIES EXPORT_NAME headers)
target_compile_features(alps_headers INTERFACE cxx_std_17)
target_compile_features(alps PRIVATE c_std_11)
target_include_directories(alps_headers INTERFACE
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/bindings>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(alps_headers INTERFACE Boost::headers)
set(ALPS_HEADER_DEPENDENCIES ${ALPS_HDF5_TARGETS})
if(ALPS_ENABLE_MPI)
  list(APPEND ALPS_HEADER_DEPENDENCIES MPI::MPI_CXX)
endif()
if(ALPS_ENABLE_OPENMP)
  list(APPEND ALPS_HEADER_DEPENDENCIES OpenMP::OpenMP_CXX)
endif()
foreach(dependency IN LISTS ALPS_HEADER_DEPENDENCIES)
  foreach(property IN ITEMS INCLUDE_DIRECTORIES COMPILE_DEFINITIONS COMPILE_OPTIONS COMPILE_FEATURES)
    set_property(TARGET alps_headers APPEND PROPERTY INTERFACE_${property}
      "$<TARGET_PROPERTY:${dependency},INTERFACE_${property}>")
  endforeach()
endforeach()
target_compile_definitions(alps PRIVATE
  BOOST_TIMER_ENABLE_DEPRECATED BOOST_ALLOW_DEPRECATED_HEADERS)
target_link_libraries(alps PUBLIC alps_headers ${Boost_LIBRARIES} ${ALPS_HDF5_TARGETS} Threads::Threads)
if(ALPS_BLAS_TARGET)
  target_link_libraries(alps PUBLIC ${ALPS_BLAS_TARGET})
endif()
if(LAPACK_FOUND)
  target_link_libraries(alps PUBLIC ${ALPS_LAPACK_TARGET})
endif()
if(BIND_FORTRAN_LOWERCASE)
  target_compile_definitions(alps_headers INTERFACE BIND_FORTRAN_LOWERCASE)
endif()
if(BIND_FORTRAN_INTEGER_8)
  target_compile_definitions(alps_headers INTERFACE BIND_FORTRAN_INTEGER_8=1)
endif()
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps PUBLIC MPI::MPI_CXX)
endif()
if(ALPS_ENABLE_OPENMP)
  target_link_libraries(alps PUBLIC OpenMP::OpenMP_CXX)
endif()
include(GenerateExportHeader)
generate_export_header(alps EXPORT_MACRO_NAME ALPS_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/src/alps/export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_STATIC_DEFINE)
endif()
set_target_properties(alps PROPERTIES
  POSITION_INDEPENDENT_CODE ON VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
install(FILES "${PROJECT_BINARY_DIR}/src/alps/export.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/alps COMPONENT headers)

if(ALPS_BUILD_FORTRAN)
  set(ALPS_FORTRAN_SOURCES fortran/fwrapper_impl.C)
  add_library(alps_fortran STATIC ${ALPS_FORTRAN_SOURCES})
  target_link_libraries(alps_fortran PUBLIC alps)
endif()
#######################################################################
# install
#######################################################################
install(TARGETS alps EXPORT ALPSTargets RUNTIME_DEPENDENCY_SET alps_runtime COMPONENT libraries
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS alps_headers EXPORT ALPSTargets COMPONENT headers)

if(ALPS_BUILD_FORTRAN)
  install(TARGETS alps_fortran EXPORT ALPSTargets COMPONENT libraries
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif(ALPS_BUILD_FORTRAN)
