#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

#
# libalps
#

add_library(alps
  utility/encode.cpp
  utility/copyright.cpp
  utility/os.cpp
  utility/temporary_filename.cpp
  utility/vmusage.cpp
  alea/nan.C
  alea/observable.C
  alea/observablefactory.C
  alea/observableset.C
  expression/evaluator.C
  lattice/disorder.C
  lattice/latticedescriptor.C
  lattice/latticegraphdescriptor.C
  lattice/latticelibrary.C
  lattice/unitcell.C
  model/bondterm.C
  model/globaloperator.C
  model/modellibrary.C
  model/siteterm.C
  osiris/buffer.C
  osiris/dump.C
  osiris/dumparchive.C
  osiris/process.C
  osiris/comm.C
  osiris/mpdump.C
  osiris/xdrdump.C
  parameter/parameter.C
  parameter/parameterlist.C
  parameter/parameters.C
  parser/parser.C
  parser/xmlattributes.C
  parser/xmlhandler.C
  parser/xmlparser.C
  parser/xmlstream.C
  parser/xslt_path.C
  hdf5/archive.cpp
  ngs/lib/signal.cpp
  ngs/lib/stacktrace.cpp
  random/rngfactory.C
  random/parallel/get_prime.cpp
  random/parallel/primelist_64.cpp
  scheduler/abstract_task.C
  scheduler/convertxdr.C
  scheduler/factory.C
  scheduler/info.C
  scheduler/master_scheduler.C
  scheduler/measurement_operators.C
  scheduler/montecarlo.C
  scheduler/mpp_scheduler.C
  scheduler/options.C
  scheduler/remote_task.C
  scheduler/remote_worker.C
  scheduler/scheduler.C
  scheduler/serial_scheduler.C
  scheduler/signal.C
  scheduler/single_scheduler.C
  scheduler/slave_task.C
  scheduler/task.C
  scheduler/worker.C
  scheduler/workertask.C
  stop_callback.cpp
  ngs/lib/api.cpp
  ngs/lib/short_print.cpp
  ngs/lib/sleep.cpp
  ngs/lib/paramvalue.cpp
  ngs/lib/paramproxy.cpp
  ngs/lib/params.cpp
  mcbase.cpp
  parseargs.cpp
  ngs/lib/mcoptions.cpp
  ngs/lib/observablewrappers.cpp
  ngs/lib/make_deprecated_parameters.cpp
  ngs/lib/make_parameters_from_xml.cpp
  ngs/lib/mcresult.cpp
  ngs/lib/mcresults.cpp
  ngs/lib/mcobservables.cpp
  ngs/lib/mcobservable.cpp
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_array.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_float.c"
  "${PROJECT_SOURCE_DIR}/third_party/xdr/src/xdr_stdio.c"
  parapack/clone.C
  parapack/clone_info.C
  parapack/filelock.C
  parapack/logger.C
  parapack/job.C
  parapack/mc_worker.C
  parapack/measurement.C
  parapack/process_impl.C
  parapack/option.C
  parapack/parapack.C
  parapack/queue.C
  parapack/rng_helper.C
  parapack/types.C
  parapack/util.C
  parapack/version.C
  parapack/worker_factory.C
  parapack/staging.C
  ngs/lib/clone.cpp
  ngs/lib/clone_info.cpp
  ngs/lib/job.cpp
  ngs/lib/parapack.cpp
  ngs/lib/worker_factory.cpp)


add_library(ALPS::alps ALIAS alps)
add_library(alps_headers INTERFACE)
add_library(ALPS::headers ALIAS alps_headers)
set_target_properties(alps_headers PROPERTIES EXPORT_NAME headers)
target_compile_features(alps_headers INTERFACE cxx_std_17)
target_compile_features(alps PRIVATE c_std_11)
target_include_directories(
  alps_headers
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/include>
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third_party/xdr/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(alps_headers INTERFACE Boost::headers)
set(ALPS_HEADER_DEPENDENCIES ${ALPS_HDF5_TARGETS})
if(ALPS_ENABLE_MPI)
  list(APPEND ALPS_HEADER_DEPENDENCIES MPI::MPI_CXX)
endif()
if(ALPS_ENABLE_OPENMP)
  list(APPEND ALPS_HEADER_DEPENDENCIES OpenMP::OpenMP_CXX)
endif()
target_link_libraries(alps_headers INTERFACE "$<COMPILE_ONLY:${ALPS_HEADER_DEPENDENCIES}>")
target_compile_definitions(alps PRIVATE BOOST_TIMER_ENABLE_DEPRECATED
                                        BOOST_ALLOW_DEPRECATED_HEADERS)
target_link_libraries(alps PUBLIC alps_headers ${ALPS_HDF5_TARGETS} Threads::Threads
                                  ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
foreach(component IN LISTS ALPS_BOOST_COMPONENTS)
  target_link_libraries(alps PUBLIC Boost::${component})
endforeach()
if(ALPS_ENABLE_MPI)
  target_link_libraries(alps PUBLIC MPI::MPI_CXX)
endif()
if(ALPS_ENABLE_OPENMP)
  target_link_libraries(alps PUBLIC OpenMP::OpenMP_CXX)
endif()
include(GenerateExportHeader)
generate_export_header(alps EXPORT_MACRO_NAME ALPS_DECL EXPORT_FILE_NAME
                       "${PROJECT_BINARY_DIR}/src/alps/export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_headers INTERFACE ALPS_STATIC_DEFINE)
endif()
set_target_properties(
  alps
  PROPERTIES VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR}
             INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_sources(alps_headers INTERFACE
  FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/src"
  FILES "${PROJECT_BINARY_DIR}/src/alps/config.h"
        "${PROJECT_BINARY_DIR}/src/alps/version.h"
        "${PROJECT_BINARY_DIR}/src/alps/export.h")

# This is a C++ bridge; only its consumers need a Fortran compiler.
add_library(alps_fortran STATIC fortran/fwrapper_impl.C)
add_library(ALPS::fortran ALIAS alps_fortran)
set_target_properties(alps_fortran PROPERTIES EXPORT_NAME fortran)
target_link_libraries(alps_fortran PUBLIC alps)
# The F77 entry points accept several argument types.
target_compile_options(alps_fortran INTERFACE
  "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-fallow-argument-mismatch>")
install(
  TARGETS alps
  EXPORT ALPSTargets
  RUNTIME_DEPENDENCY_SET alps_runtime
  COMPONENT libraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  TARGETS alps_headers alps_fortran
  EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
  FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/alps"
            "${PROJECT_SOURCE_DIR}/src/ietl"
            "${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/include/boost"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT headers
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.ipp"
  PATTERN "resources" EXCLUDE)

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/third_party/xdr/include/"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT headers)
install(
  FILES "${PROJECT_SOURCE_DIR}/third_party/xdr/LICENSE"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/xdr
  COMPONENT libraries)

install(
  FILES "${PROJECT_SOURCE_DIR}/third_party/boost_numeric_bindings/LICENSE_1_0.txt"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/licenses/boost_numeric_bindings
  COMPONENT headers)

install(
  DIRECTORY resources/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/xml
  COMPONENT xml
  FILES_MATCHING PATTERN "*.xml" PATTERN "*.xsl")
