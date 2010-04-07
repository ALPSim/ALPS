##########################################################################
# Core Functionality for Boost                                           #
##########################################################################
# Copyright (C) 2007-2008 Douglas Gregor <doug.gregor@gmail.com>         #
# Copyright (C) 2007 Troy Straszheim                                     #
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.             #
# See accompanying file LICENSE_1_0.txt or copy at                       #
#   http://www.boost.org/LICENSE_1_0.txt                                 #
##########################################################################
# Important developer macros in this file:                               #
#                                                                        #
#   boost_library_project: Defines a Boost library project (e.g.,        #
#   Boost.Python).                                                       #
#                                                                        #
#   boost_add_library: Builds library binaries for Boost libraries       #
#   with compiled sources (e.g., boost_filesystem).                      #
#                                                                        #
#   boost_add_executable: Builds executables.                            #
##########################################################################

# Creates a new executable from source files.
#
#   boost_add_executable(exename
#                        source1 source2 ...
#                        [COMPILE_FLAGS compileflags]
#                        [feature_COMPILE_FLAGS compileflags]
#                        [LINK_FLAGS linkflags]
#                        [feature_LINK_FLAGS linkflags]
#                        [LINK_LIBS linklibs]
#                        [feature_LINK_LIBS linklibs]
#                        [DEPENDS libdepend1 libdepend2 ...]
#                        [feature]
#                        [NO_INSTALL])
#
# where exename is the name of the executable (e.g., "wave").  source1,
# source2, etc. are the source files used to build the executable, e.g.,
# cpp.cpp. If no source files are provided, "exename.cpp" will be
# used.
#
# This macro has a variety of options that affect its behavior. In
# several cases, we use the placeholder "feature" in the option name
# to indicate that there are actually several different kinds of
# options, each referring to a different build feature, e.g., shared
# libraries, multi-threaded, debug build, etc. For a complete listing
# of these features, please refer to the CMakeLists.txt file in the
# root of the Boost distribution, which defines the set of features
# that will be used to build Boost libraries by default.
#
# The options that affect this macro's behavior are:
#
#   COMPILE_FLAGS: Provides additional compilation flags that will be
#   used when building the executable.
#
#   feature_COMPILE_FLAGS: Provides additional compilation flags that
#   will be used only when building the executable with the given
#   feature (e.g., SHARED_COMPILE_FLAGS when we're linking against
#   shared libraries). Note that the set of features used to build the
#   executable depends both on the arguments given to
#   boost_add_executable (see the "feature" argument description,
#   below) and on the user's choice of variants to build.
#
#   LINK_FLAGS: Provides additional flags that will be passed to the
#   linker when linking the executable. This option should not be used
#   to link in additional libraries; see LINK_LIBS and DEPENDS.
#
#   feature_LINK_FLAGS: Provides additional flags that will be passed
#   to the linker when linking the executable with the given feature
#   (e.g., MULTI_THREADED_LINK_FLAGS when we're linking a
#   multi-threaded executable).
#
#   LINK_LIBS: Provides additional libraries against which the
#   executable will be linked. For example, one might provide "expat"
#   as options to LINK_LIBS, to state that the executable will link
#   against the expat library binary. Use LINK_LIBS for libraries
#   external to Boost; for Boost libraries, use DEPENDS.
#
#   feature_LINK_LIBS: Provides additional libraries to link against
#   when linking an executable built with the given feature. 
#
#   DEPENDS: States that this executable depends on and links against
#   a Boostlibrary. The arguments to DEPENDS should be the unversioned
#   name of the Boost library, such as "boost_filesystem". Like
#   LINK_LIBS, this option states that the executable will link
#   against the stated libraries. Unlike LINK_LIBS, however, DEPENDS
#   takes particular library variants into account, always linking to
#   the appropriate variant of a Boost library. For example, if the
#   MULTI_THREADED feature was requested in the call to
#   boost_add_executable, DEPENDS will ensure that we only link
#   against multi-threaded libraries.
#
#   feature: States that the executable should always be built using a
#   given feature, e.g., SHARED linking (against its libraries) or
#   MULTI_THREADED (for multi-threaded builds). If that feature has
#   been turned off by the user, the executable will not build.
#
#   NO_INSTALL: Don't install this executable with the rest of Boost.
#
#   OUTPUT_NAME: If you want the executable to be generated somewhere
#   other than the binary directory, pass the path (including
#   directory and file name) via the OUTPUT_NAME parameter.
#
# Example:
#   boost_add_executable(wave cpp.cpp 
#     DEPENDS boost_wave boost_program_options boost_filesystem 
#             boost_serialization
#     )
macro(boost_add_executable EXENAME)
  # Note: ARGS is here to support the use of boost_add_executable in
  # the testing code.
  parse_arguments(THIS_EXE
    "DEPENDS;COMPILE_FLAGS;LINK_FLAGS;LINK_LIBS;OUTPUT_NAME;ARGS;${BOOST_ADD_ARG_NAMES}"
    "NO_INSTALL;${BOOST_ADDEXE_OPTION_NAMES}"
    ${ARGN}
    )

  # Determine the list of sources
  if (THIS_EXE_DEFAULT_ARGS)
    set(THIS_EXE_SOURCES ${THIS_EXE_DEFAULT_ARGS})
  else (THIS_EXE_DEFAULT_ARGS)
    set(THIS_EXE_SOURCES ${EXENAME}.cpp)
  endif (THIS_EXE_DEFAULT_ARGS)
  

    # Compute the actual set of library dependencies, based on the
    # variant name we computed above. The RELEASE and DEBUG versions
    # only apply when THIS_EXE_DEBUG_AND_RELEASE.
    set(THIS_EXE_ACTUAL_DEPENDS)
    foreach(LIB ${THIS_EXE_DEPENDS})
        list(APPEND THIS_EXE_ACTUAL_DEPENDS "${LIB}${VARIANT_TARGET_NAME}")
    endforeach(LIB ${THIS_EXE_DEPENDS})

    add_executable(${EXENAME} ${THIS_EXE_SOURCES})
    
    # Set the various compilation and linking flags
    set_target_properties(${EXENAME}
      PROPERTIES
      COMPILE_FLAGS "${THIS_EXE_COMPILE_FLAGS}"
      LINK_FLAGS "${THIS_EXE_LINK_FLAGS}"
      )

    # If the user gave an output name, use it.
    if(THIS_EXE_OUTPUT_NAME)
      set_target_properties(${EXENAME}
        PROPERTIES
        OUTPUT_NAME ${THIS_EXE_OUTPUT_NAME}
        )
    endif()

    # Link against the various libraries 
      target_link_libraries(${EXENAME} 
        ${THIS_EXE_ACTUAL_DEPENDS} 
        ${THIS_EXE_LINK_LIBS})

    # Install the executable, if not suppressed
    if (NOT THIS_EXE_NO_INSTALL)
      install(TARGETS ${EXENAME} DESTINATION bin)
    endif (NOT THIS_EXE_NO_INSTALL)
endmacro(boost_add_executable)

macro(alps_library_project LIBNAME)
  project(${LIBNAME})
  if(EXISTS ${alps2_SOURCE_DIR}/libs/${PROJECT_NAME}/src)
      add_subdirectory(src)
  endif()
  
  if(BUILD_TESTING)
	if(EXISTS  ${alps2_SOURCE_DIR}/libs/${PROJECT_NAME}/test)
	  # Testing is enabled globally and this project has some
	  # tests. So, include the tests
	  add_custom_target(${PROJECT_NAME}-test)
	  add_dependencies(test ${PROJECT_NAME}-test)
      add_subdirectory(test)
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/tests)
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)
      endif(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/tests)
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/tests/${PROJECT_NAME})
        file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/${PROJECT_NAME})
      endif(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/tests/${PROJECT_NAME})
	endif()
  endif(BUILD_TESTING)

  if (BUILD_DOCUMENTATION)
    if(EXISTS  ${alps2_SOURCE_DIR}/libs/${PROJECT_NAME}/doc)
      add_subdirectory(doc)
    endif ()
  endif()

endmacro(alps_library_project)





