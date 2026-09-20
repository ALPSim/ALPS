# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.22)
# The installed SDK supplies one runtime directory. Per-module vcpkg copies
# would introduce duplicate DLL paths when CMake collects wheel dependencies.
set(VCPKG_APPLOCAL_DEPS OFF)
project(pyalps LANGUAGES C CXX)

option(PYALPS_BUILD_APPLICATIONS "Build optional ALPS application bindings" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ALPS REQUIRED CONFIG)
if(NOT ALPS_BUILD_SHARED_LIBS)
  message(FATAL_ERROR "pyalps requires a shared ALPS SDK so every extension uses the same C++ runtime")
endif()
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

# Policy (see README "Free-threading and stable-ABI policy"): do NOT
# add FREE_THREADED (libalps relies on the GIL as its lock around
# shared state) and do not add STABLE_ABI without revisiting the wheel
# matrix — per-version wheels are deliberate.
# NB_STATIC keeps nanobind's core private to each extension. On musllinux,
# MUSL_DYNAMIC_LIBCPP prevents nanobind from also embedding libstdc++ and
# libgcc into every module; auditwheel vendors those shared runtimes once.
set(_pyalps_nanobind_options NB_STATIC MUSL_DYNAMIC_LIBCPP)

nanobind_add_module(pyalea_c ${_pyalps_nanobind_options} "${_bindings}/pyalea.cpp")
nanobind_add_module(pymcdata_c ${_pyalps_nanobind_options} "${_bindings}/pymcdata.cpp")
nanobind_add_module(pytools_c ${_pyalps_nanobind_options} "${_bindings}/pytools.cpp")
nanobind_add_module(pyngsparams_c ${_pyalps_nanobind_options} "${_bindings}/ngs/params.cpp")
nanobind_add_module(pyngshdf5_c ${_pyalps_nanobind_options} "${_bindings}/ngs/hdf5.cpp")
nanobind_add_module(pyngsbase_c ${_pyalps_nanobind_options} "${_bindings}/ngs/mcbase.cpp")
nanobind_add_module(pyngsobservable_c ${_pyalps_nanobind_options} "${_bindings}/ngs/observable.cpp")
nanobind_add_module(pyngsobservables_c ${_pyalps_nanobind_options} "${_bindings}/ngs/observables.cpp")
nanobind_add_module(pyngsresult_c ${_pyalps_nanobind_options} "${_bindings}/ngs/result.cpp")
nanobind_add_module(pyngsresults_c ${_pyalps_nanobind_options} "${_bindings}/ngs/results.cpp")
nanobind_add_module(pyngsapi_c ${_pyalps_nanobind_options} "${_bindings}/ngs/api.cpp")
nanobind_add_module(pyngsrandom01_c ${_pyalps_nanobind_options} "${_bindings}/ngs/random01.cpp")
nanobind_add_module(pyngsaccumulator_c ${_pyalps_nanobind_options} "${_bindings}/ngs/accumulator.cpp")

if(PYALPS_BUILD_APPLICATIONS)
  if(NOT EXISTS "${_alps_source_root}/tool/maxent.cpp")
    message(FATAL_ERROR
      "PYALPS_BUILD_APPLICATIONS requires an ALPS source checkout or pyalps sdist. "
      "Configure with -DPYALPS_BUILD_APPLICATIONS=OFF for the core-only package.")
  endif()

  set(_dmft "${_alps_source_root}/applications/dmft/qmc")

  nanobind_add_module(maxent_c ${_pyalps_nanobind_options}
    "${_alps_source_root}/tool/maxent.cpp"
    "${_alps_source_root}/tool/maxent_helper.cpp"
    "${_alps_source_root}/tool/maxent_simulation.cpp"
    "${_alps_source_root}/tool/maxent_parms.cpp")

  nanobind_add_module(cthyb ${_pyalps_nanobind_options}
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

  nanobind_add_module(ctint ${_pyalps_nanobind_options}
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

  list(APPEND _pyalps_targets maxent_c cthyb ctint)
  foreach(_target IN ITEMS maxent_c cthyb ctint)
    target_compile_definitions(${_target} PRIVATE BUILD_PYTHON_MODULE)
    target_include_directories(${_target} PRIVATE "${_bindings}" "${_dmft}")
  endforeach()
endif()

option(PYALPS_BUNDLE_APPLICATIONS
  "Bundle the ALPS application executables into pyalps/bin" ON)

# Where the extension modules look for libalps at run time.
#
# When the ALPS libraries ship inside the package (the default), the *only*
# route to libalps must be the in-package pyalps/lib copy -- the same one the
# bundled programs load. Leaving the SDK's own library directory on the list
# as well is not a harmless fallback: auditwheel/delocate resolve the first
# match, find libalps outside the wheel, and vendor a second copy into
# pyalps.libs (or pyalps/.dylibs), so the wheel carries the same 10 MB library
# twice. HDF5 stays on the list precisely because it *should* be vendored.
#
# With PYALPS_BUNDLE_APPLICATIONS=OFF nothing installs libalps into the
# package, so there the SDK directory is the only way for the repair tools to
# find it at all.
if(APPLE)
  set(_pyalps_in_package_libs "@loader_path/../lib")
else()
  set(_pyalps_in_package_libs "$ORIGIN/../lib")
endif()
set(_pyalps_hdf5_runtime_dirs "")
foreach(_library IN LISTS HDF5_LIBRARIES)
  if(TARGET "${_library}")
    get_target_property(_type "${_library}" TYPE)
    if(_type MATCHES "^(SHARED|UNKNOWN)_LIBRARY$")
      list(APPEND _pyalps_hdf5_runtime_dirs "$<TARGET_FILE_DIR:${_library}>")
    endif()
  elseif(IS_ABSOLUTE "${_library}")
    get_filename_component(_directory "${_library}" DIRECTORY)
    list(APPEND _pyalps_hdf5_runtime_dirs "${_directory}")
  endif()
endforeach()
list(REMOVE_DUPLICATES _pyalps_hdf5_runtime_dirs)
if(PYALPS_BUNDLE_APPLICATIONS)
  set(_pyalps_install_rpath
    "${_pyalps_in_package_libs};${_pyalps_hdf5_runtime_dirs}")
else()
  set(_pyalps_install_rpath
    "${ALPS_LIBRARY_DIRS};${_pyalps_hdf5_runtime_dirs}")
endif()

foreach(_target IN LISTS _pyalps_targets)
  target_link_libraries(${_target} PRIVATE ALPS::alps)
  set_target_properties(${_target} PROPERTIES
    INSTALL_RPATH "${_pyalps_install_rpath}")
endforeach()

set(_extension_dir "pyalps/_ext")
install(TARGETS ${_pyalps_targets}
  RUNTIME_DEPENDENCY_SET pyalps_runtime
  LIBRARY DESTINATION "${_extension_dir}"
  RUNTIME DESTINATION "${_extension_dir}")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/pyalps/"
  DESTINATION pyalps
  FILES_MATCHING PATTERN "*.py")

# Runtime fallback paths for pyalps.tools, used only when the package does
# not carry the resource in question. The XML/XSL library below is always
# bundled, so its fallback is always empty; the application binaries are
# bundled unless PYALPS_BUNDLE_APPLICATIONS is OFF, and only that
# configuration needs to point back at the SDK that built the package.
#
# Never fill these in for a bundled package: they would be absolute paths on
# the build machine, baked into a redistributable artifact.
set(PYALPS_ALPS_XML_FALLBACK "")
if(PYALPS_BUNDLE_APPLICATIONS)
  set(PYALPS_ALPS_BIN_FALLBACK "")
else()
  set(PYALPS_ALPS_BIN_FALLBACK "${ALPS_ROOT_DIR}/bin")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/pyalps/pyalps_config.py.in"
  "${CMAKE_CURRENT_BINARY_DIR}/pyalps_config.py" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pyalps_config.py" DESTINATION pyalps)

# The ALPS application executables ship inside the package (pyalps/bin), as
# the legacy ALPS_PYTHON_WHEEL build did. Keep this list explicit: the SDK's
# bin directory also contains conversion and plotting tools whose relative
# data paths assume a prefix install rather than the wheel layout.
# pyalps.tools.check_existence() still prepends <pyalps>/bin to PATH, so
# without these files runApplication('spinmc', ...) fails on a wheel install
# with "There is no spinmc on the path!" -- the fallback it tries instead,
# pyalps_config.ALPS_BIN_INSTALL_DIR, is the build machine's SDK path.
#
# The programs are copied from the SDK rather than rebuilt here. The root
# CMakeLists gives every installed target an RPATH of @loader_path/../lib
# (APPLE) or $ORIGIN/../lib, so they resolve libalps out of pyalps/lib below
# without any post-install patching -- the same directory the extension
# modules point at, so one copy of each library serves both.
if(PYALPS_BUNDLE_APPLICATIONS)
  set(_alps_bin_dir "${ALPS_ROOT_DIR}/bin")
  set(_alps_application_names
    checksign
    dirloop_sse
    dmft
    dmrg
    fulldiag
    fulldiag_evaluate
    hirschfye
    hybridization
    interaction
    loop
    qwl
    qwl_evaluate
    simplemc
    sparsediag
    spinmc
    spinmc_evaluate
    worm
    worm_evaluate)
  set(_alps_programs "")
  set(_alps_program_targets "")
  set(_alps_missing_programs "")
  foreach(_program IN LISTS _alps_application_names)
    if(EXISTS "${_alps_bin_dir}/${_program}${CMAKE_EXECUTABLE_SUFFIX}")
      list(APPEND _alps_programs "${_alps_bin_dir}/${_program}${CMAKE_EXECUTABLE_SUFFIX}")
      add_executable(_pyalps_program_${_program} IMPORTED)
      set_target_properties(_pyalps_program_${_program} PROPERTIES
        IMPORTED_LOCATION "${_alps_bin_dir}/${_program}${CMAKE_EXECUTABLE_SUFFIX}")
      list(APPEND _alps_program_targets _pyalps_program_${_program})
    else()
      list(APPEND _alps_missing_programs "${_program}")
    endif()
  endforeach()
  if(_alps_missing_programs)
    message(FATAL_ERROR
      "PYALPS_BUNDLE_APPLICATIONS is ON but the wheel SDK is missing these "
      "legacy wheel applications: ${_alps_missing_programs}. Build the ALPS "
      "SDK with the wheel-deps preset, or configure with "
      "-DPYALPS_BUNDLE_APPLICATIONS=OFF for a bindings-only wheel -- in which "
      "case pyalps.runApplication() will not work.")
  endif()
  list(LENGTH _alps_programs _alps_program_count)
  message(STATUS "pyalps: bundling ${_alps_program_count} ALPS programs into pyalps/bin")
  install(IMPORTED_RUNTIME_ARTIFACTS ${_alps_program_targets}
    RUNTIME_DEPENDENCY_SET pyalps_runtime RUNTIME DESTINATION pyalps/bin)

  # Imported target metadata identifies the actual runtime and its SONAME on
  # ELF and macOS alike. Copy its bytes once, using the loader-visible name.
  if(NOT WIN32)
    install(FILES "$<TARGET_FILE:ALPS::alps>" DESTINATION pyalps/lib
      RENAME "$<TARGET_SONAME_FILE_NAME:ALPS::alps>")
  endif()
endif()

# Windows extensions and child processes share the DLLs next to the bundled
# executables. Python registers this directory before loading any extension.
if(WIN32)
  install(RUNTIME_DEPENDENCY_SET pyalps_runtime
    DIRECTORIES "${ALPS_RUNTIME_LIBRARY_DIRS}"
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*" "python[0-9]+.*\\.dll"
    POST_EXCLUDE_REGEXES ".*[/\\\\][Ww][Ii][Nn][Dd][Oo][Ww][Ss][/\\\\].*"
    RUNTIME DESTINATION pyalps/bin)
endif()
include("${_alps_source_root}/cmake/ALPSRuntime.cmake")
alps_install_compiler_runtime(pyalps/bin)

# Run after all runtime files have been installed. Wheel repair regenerates this
# manifest once auditwheel/delocate have finalized their paths and install names.
install(CODE "execute_process(COMMAND \"${Python_EXECUTABLE}\"
  \"${CMAKE_CURRENT_SOURCE_DIR}/_build_support/runtime_manifest.py\"
  --tree \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}\" --alps-version \"${ALPS_VERSION}\"
  COMMAND_ERROR_IS_FATAL ANY)")

# The ALPS XML/XSL library ships inside the package (pyalps/xml), as the
# legacy wheel build did: the stylesheets plus the lattice and model
# libraries that parameter files reference by default.
set(_alps_xml_source "${_alps_source_root}/lib/xml")
if(NOT EXISTS "${_alps_xml_source}/ALPS.xsl")
  message(FATAL_ERROR
    "ALPS XML stylesheets not found at ${_alps_xml_source}. "
    "They are required so the pyalps package can bundle pyalps/xml.")
endif()
install(DIRECTORY "${_alps_xml_source}/" DESTINATION pyalps/xml
  FILES_MATCHING PATTERN "*.xsl")
configure_file("${_alps_xml_source}/lattices.xml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xml/lattices.xml" COPYONLY)
configure_file("${_alps_xml_source}/models.xml.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xml/models.xml" COPYONLY)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/xml/lattices.xml"
  "${CMAKE_CURRENT_BINARY_DIR}/xml/models.xml"
  DESTINATION pyalps/xml)
