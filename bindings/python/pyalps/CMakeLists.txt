# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.22)
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

# Compile the bindings with the same preprocessor configuration the
# ALPS SDK was built with (BOOST_NO_AUTO_PTR and friends), so the Boost
# headers on both sides of the library boundary are configured
# identically. The -D entries are lifted from the exported
# ALPS_CMAKE_CXX_FLAGS rather than hand-copied from the root
# CMakeLists, so a define added there cannot silently desync.
separate_arguments(_alps_sdk_cxx_flags NATIVE_COMMAND "${ALPS_CMAKE_CXX_FLAGS}")
set(_alps_sdk_definitions "")
foreach(_flag IN LISTS _alps_sdk_cxx_flags)
  if(_flag MATCHES "^-D(.+)$")
    list(APPEND _alps_sdk_definitions "${CMAKE_MATCH_1}")
  endif()
endforeach()

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
if(PYALPS_BUNDLE_APPLICATIONS)
  set(_pyalps_install_rpath
    "${_pyalps_in_package_libs};${_alps_hdf5_prefix}/lib")
else()
  set(_pyalps_install_rpath
    "${ALPS_LIBRARY_DIRS};${_alps_hdf5_prefix}/lib")
endif()

foreach(_target IN LISTS _pyalps_targets)
  if(_alps_sdk_definitions)
    target_compile_definitions(${_target} PRIVATE ${_alps_sdk_definitions})
  endif()
  target_include_directories(${_target} PRIVATE
    ${ALPS_INCLUDE_DIRS}
    ${ALPS_EXTRA_INCLUDE_DIRS})
  target_link_libraries(${_target} PRIVATE ${_pyalps_link_libraries})
  target_link_options(${_target} PRIVATE ${_pyalps_link_options})
  set_target_properties(${_target} PROPERTIES
    INSTALL_RPATH "${_pyalps_install_rpath}")
endforeach()

set(_extension_dir "pyalps/_ext")
install(TARGETS ${_pyalps_targets}
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
  set(_alps_missing_programs "")
  foreach(_program IN LISTS _alps_application_names)
    if(EXISTS "${_alps_bin_dir}/${_program}")
      list(APPEND _alps_programs "${_alps_bin_dir}/${_program}")
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
  install(PROGRAMS ${_alps_programs} DESTINATION pyalps/bin)

  # The bundled programs resolve libalps and libboost through their ../lib
  # RPATH, so the SDK's shared libraries have to travel with them -- each one
  # exactly once. install(DIRECTORY) is the wrong tool here: it would copy the
  # whole SONAME chain (libalps.so, libalps.so.2, libalps.so.2.3.4), and since
  # a wheel cannot carry symlinks those arrive as three byte-identical 10 MB
  # files. Install only the SONAME-named file, which is the name every
  # DT_NEEDED / LC_LOAD_DYLIB entry actually asks for.
  #
  # src/alps/CMakeLists.txt gives libalps SOVERSION ${ALPS_VERSION_MAJOR} on
  # ELF platforms and no version at all on APPLE, so "<name>.<major>" is the
  # SONAME spelling wherever a versioned form exists. test/pyalps guards both
  # halves of this: that no library is bundled twice, and that every bundled
  # program can still be loaded.
  file(GLOB _alps_runtime_libs LIST_DIRECTORIES FALSE
    "${ALPS_ROOT_DIR}/lib/*.so" "${ALPS_ROOT_DIR}/lib/*.dylib")
  if(NOT _alps_runtime_libs)
    message(FATAL_ERROR
      "PYALPS_BUNDLE_APPLICATIONS is ON but no shared libraries were found in "
      "${ALPS_ROOT_DIR}/lib. The bundled programs would ship without the "
      "libalps they link against.")
  endif()
  foreach(_lib IN LISTS _alps_runtime_libs)
    get_filename_component(_lib_name "${_lib}" NAME)
    set(_lib_soname "${_lib_name}")
    if(EXISTS "${_lib}.${ALPS_VERSION_MAJOR}")
      set(_lib_soname "${_lib_name}.${ALPS_VERSION_MAJOR}")
    endif()
    get_filename_component(_lib_real "${_lib}" REALPATH)
    install(PROGRAMS "${_lib_real}" DESTINATION pyalps/lib
      RENAME "${_lib_soname}")
  endforeach()
endif()

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
