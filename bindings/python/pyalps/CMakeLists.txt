# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.18)
project(pyalps LANGUAGES CXX)

option(PYALPS_BUILD_APPLICATIONS "Build optional ALPS application bindings" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ALPS REQUIRED CONFIG)
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _nanobind_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
list(PREPEND CMAKE_PREFIX_PATH "${_nanobind_cmake_dir}")
find_package(nanobind 2.10 CONFIG REQUIRED)

get_filename_component(_repo_root "${CMAKE_CURRENT_SOURCE_DIR}/../../.." ABSOLUTE)
set(_alps_source_root "${_repo_root}")
if(NOT EXISTS "${_alps_source_root}/tool/maxent.cpp")
  set(_alps_source_root "${CMAKE_CURRENT_SOURCE_DIR}/_vendor")
endif()
set(_bindings "${CMAKE_CURRENT_SOURCE_DIR}/cpp")

link_directories(${ALPS_LIBRARY_DIRS})
if(ALPS_HDF5_INCLUDE_DIR)
  get_filename_component(_alps_hdf5_prefix "${ALPS_HDF5_INCLUDE_DIR}" DIRECTORY)
  link_directories("${_alps_hdf5_prefix}/lib")
endif()

# Legacy ALPSConfig.cmake serializes the HDF5 imported target name into
# ALPS_LIBRARIES.  Outside the original build tree that target no longer
# exists, while the installed library keeps its conventional `hdf5` name.
set(_pyalps_link_libraries ${ALPS_LIBRARIES})
list(TRANSFORM _pyalps_link_libraries REPLACE "^hdf5-shared$" "hdf5")
separate_arguments(_pyalps_link_options NATIVE_COMMAND "${ALPS_EXTRA_LINKER_FLAGS}")

set(_pyalps_targets
  pyalea_c
  pymcdata_c
  pytools_c
  pyngsparams_c
  pyngshdf5_c
  pyngsbase_c
  pyngsobservable_c
  pyngsobservables_c
  pyngsresult_c
  pyngsresults_c
  pyngsapi_c
  pyngsrandom01_c
  pyngsaccumulator_c)

nanobind_add_module(pyalea_c NB_STATIC "${_bindings}/pyalea.cpp")
nanobind_add_module(pymcdata_c NB_STATIC "${_bindings}/pymcdata.cpp")
nanobind_add_module(pytools_c NB_STATIC "${_bindings}/pytools.cpp")
nanobind_add_module(pyngsparams_c NB_STATIC "${_bindings}/ngs/params.cpp")
nanobind_add_module(pyngshdf5_c NB_STATIC "${_bindings}/ngs/hdf5.cpp")
nanobind_add_module(pyngsbase_c NB_STATIC "${_bindings}/ngs/mcbase.cpp")
nanobind_add_module(pyngsobservable_c NB_STATIC "${_bindings}/ngs/observable.cpp")
nanobind_add_module(pyngsobservables_c NB_STATIC "${_bindings}/ngs/observables.cpp")
nanobind_add_module(pyngsresult_c NB_STATIC "${_bindings}/ngs/result.cpp")
nanobind_add_module(pyngsresults_c NB_STATIC "${_bindings}/ngs/results.cpp")
nanobind_add_module(pyngsapi_c NB_STATIC "${_bindings}/ngs/api.cpp")
nanobind_add_module(pyngsrandom01_c NB_STATIC "${_bindings}/ngs/random01.cpp")
nanobind_add_module(pyngsaccumulator_c NB_STATIC "${_bindings}/ngs/accumulator.cpp")

if(PYALPS_BUILD_APPLICATIONS)
  if(NOT EXISTS "${_alps_source_root}/tool/maxent.cpp")
    message(FATAL_ERROR
      "PYALPS_BUILD_APPLICATIONS requires an ALPS source checkout or pyalps sdist. "
      "Configure with -DPYALPS_BUILD_APPLICATIONS=OFF for the core-only package.")
  endif()

  set(_dmft "${_alps_source_root}/applications/dmft/qmc")

  nanobind_add_module(maxent_c NB_STATIC
    "${_alps_source_root}/tool/maxent.cpp"
    "${_alps_source_root}/tool/maxent_helper.cpp"
    "${_alps_source_root}/tool/maxent_simulation.cpp"
    "${_alps_source_root}/tool/maxent_parms.cpp")

  nanobind_add_module(cthyb NB_STATIC
    "${_dmft}/hybridization/hybmain.cpp"
    "${_dmft}/hybridization/hybsim.cpp"
    "${_dmft}/hybridization/hyblocal.cpp"
    "${_dmft}/hybridization/hybint.cpp"
    "${_dmft}/hybridization/hybfun.cpp"
    "${_dmft}/hybridization/hybretintfun.cpp"
    "${_dmft}/hybridization/hybmatrix.cpp"
    "${_dmft}/hybridization/hybmatrix_ft.cpp"
    "${_dmft}/hybridization/hybconfig.cpp"
    "${_dmft}/hybridization/hybupdates.cpp"
    "${_dmft}/hybridization/hybevaluate.cpp"
    "${_dmft}/hybridization/hybmeasurements.cpp")

  nanobind_add_module(ctint NB_STATIC
    "${_dmft}/interaction_expansion2/main.cpp"
    "${_dmft}/fouriertransform.C"
    "${_dmft}/interaction_expansion2/auxiliary.cpp"
    "${_dmft}/interaction_expansion2/observables.cpp"
    "${_dmft}/interaction_expansion2/fastupdate.cpp"
    "${_dmft}/interaction_expansion2/selfenergy.cpp"
    "${_dmft}/interaction_expansion2/solver.cpp"
    "${_dmft}/interaction_expansion2/io.cpp"
    "${_dmft}/interaction_expansion2/splines.cpp"
    "${_dmft}/interaction_expansion2/interaction_expansion.cpp"
    "${_dmft}/interaction_expansion2/measurements.cpp"
    "${_dmft}/interaction_expansion2/model.cpp")

  nanobind_add_module(dwa_c NB_STATIC "${_bindings}/apps/dwa.cpp")

  list(APPEND _pyalps_targets maxent_c cthyb ctint dwa_c)
  foreach(_target IN ITEMS maxent_c cthyb ctint)
    target_compile_definitions(${_target} PRIVATE BUILD_PYTHON_MODULE)
    target_include_directories(${_target} PRIVATE "${_bindings}" "${_dmft}")
  endforeach()
  target_include_directories(dwa_c PRIVATE "${_alps_source_root}/applications/qmc/dwa")
endif()

foreach(_target IN LISTS _pyalps_targets)
  target_include_directories(${_target} PRIVATE
    ${ALPS_INCLUDE_DIRS}
    ${ALPS_EXTRA_INCLUDE_DIRS})
  target_link_libraries(${_target} PRIVATE ${_pyalps_link_libraries})
  target_link_options(${_target} PRIVATE ${_pyalps_link_options})
  set_target_properties(${_target} PROPERTIES
    INSTALL_RPATH "${ALPS_LIBRARY_DIRS};${_alps_hdf5_prefix}/lib")
endforeach()

set(_extension_dir "pyalps/_ext")
install(TARGETS ${_pyalps_targets}
  LIBRARY DESTINATION "${_extension_dir}"
  RUNTIME DESTINATION "${_extension_dir}")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/pyalps/"
  DESTINATION pyalps
  FILES_MATCHING PATTERN "*.py")
