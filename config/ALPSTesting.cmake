##########################################################################
# Regression Testing Support for Boost                                   #
##########################################################################
# Copyright (C) 2007-8 Douglas Gregor <doug.gregor@gmail.com>            #
# Copyright (C) 2007-8 Troy D. Straszheim                                #
#                                                                        #
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
##########################################################################
# This file provides a set of CMake macros that support regression
# testing for Boost libraries. For each of the test macros below, the
# first argument, testname, states the name of the test that will be
# created. If no other arguments are provided, the source file
# testname.cpp will be used as the source file; otherwise, source
# files should be listed immediately after the name of the test.
#
# The macros for creating regression tests are:
#   boost_test_run: Builds an executable and runs it as a test. The test
#                   succeeds if it builds and returns 0 when executed.
#
#   boost_test_run_fail: Builds an executable and runs it as a test. The
#                        test succeeds if it builds but returns a non-zero
#                        exit code when executed.
#  
#   boost_test_compile: Tests that the given source file compiles without 
#                       any errors.
#
#   boost_test_compile_fail: Tests that the given source file produces 
#                            errors when compiled.

# User-controlled option that can be used to enable/disable regression
# testing. By default, we ena testing, because most users building from source will
# want to check whether the buils is correct
option(BUILD_TESTING "Enable testing" ON)

if (BUILD_TESTING)
  add_custom_target(test COMMENT "Running all tests")
    add_dependencies(test
      ${PROJECT_NAME}-test
      )
endif (BUILD_TESTING)


  file(TO_NATIVE_PATH ${alps2_SOURCE_DIR}/tools/build/CMake/passthru.py BOOST_TEST_DRIVER)



# This macro is an internal utility macro that helps parse the
# arguments passed to the Boost testing commands. It will generally
# not be used by Boost developers.
macro(boost_test_add_dependent_includes includes)
  foreach (include ${includes})
    #message(STATUS "include: ${include}")
    include_directories("${Boost_SOURCE_DIR}/libs/${include}/include")
  endforeach (include ${includes})
endmacro(boost_test_add_dependent_includes includes)



# This macro is an internal utility macro that helps parse the
# arguments passed to the Boost testing commands. It will generally
# not be used by Boost developers.
#
#   boost_test_parse_args(testname 
#                         [source1 source2 ...]
#                         [ARGS arg1 arg2... ]
#                         [COMPILE_FLAGS compileflags]
#                         [LINK_FLAGS linkflags]
#                         [LINK_LIBS linklibs]
#                         [DEPENDS libdepend1 libdepend2 ...]
#                         [COMPILE] [RUN] [FAIL])
#
# testname is the name of the test. The remaining arguments passed to
# this macro will be parsed and categorized for the developer-level
# test macros to use. 
#
# Variables affected:
#
#   BOOST_TEST_OKAY: Will be set to TRUE if it is okay to build and
#   run this test.
#
#   BOOST_TEST_SOURCES: Will be populated with the set of source files
#   that should be used to compile this test. If the user has provided
#   source files, BOOST_TEST_SOURCES will contain those; otherwise,
#   BOOST_TEST_SOURCES will only contain "testname.cpp".
#
#   BOOST_TEST_TAG:  compile, compile_fail, run, or run_fail.
#   Used in test-reporting systems.
#
#   BOOST_TEST_TESTNAME: A (hopefully) globally unique target name
#   for the test, constructed from PROJECT-testname-TAG
#
#   BOOST_TEST_arg: Will be populated with the arguments provided for
#   the arguemnt "arg", where "arg" can be any of the extra arguments
#   specified above.
#
#   
macro(boost_test_parse_args testname)
  #message("boost_test_parse_args ${testname} ${ARGN}")
  set(BOOST_TEST_OKAY TRUE)
  set(BOOST_TEST_COMPILE_FLAGS "")
  parse_arguments(BOOST_TEST 
    "LINK_LIBS;LINK_FLAGS;DEPENDS;COMPILE_FLAGS;ARGS;EXTRA_OPTIONS"
    "COMPILE;RUN;LINK;FAIL;RELEASE;DEBUG"
    ${ARGN}
    )
    
  # Check each of the dependencies to see if we can still build this
  # test.
  foreach(ARG ${BOOST_TEST_DEPENDS})
    get_target_property(DEPEND_TYPE ${ARG} TYPE)
    get_target_property(DEPEND_LOCATION ${ARG} LOCATION)
    # If building static libraries is turned off, don't try to build
    # the test
    if (NOT BUILD_STATIC AND ${DEPEND_TYPE} STREQUAL "STATIC_LIBRARY")
      set(BOOST_TEST_OKAY FALSE)
    endif (NOT BUILD_STATIC AND ${DEPEND_TYPE} STREQUAL "STATIC_LIBRARY")

    # If building shared libraries is turned off, don't try to build
    # the test
    if (NOT BUILD_SHARED AND ${DEPEND_TYPE} STREQUAL "SHARED_LIBRARY")
      set(BOOST_TEST_OKAY FALSE)
    endif (NOT BUILD_SHARED AND ${DEPEND_TYPE} STREQUAL "SHARED_LIBRARY")
  endforeach(ARG ${BOOST_TEST_DEPENDS})

  # Setup the SOURCES variables. If no sources are specified, use the
  # name of the test.cpp
  if (BOOST_TEST_DEFAULT_ARGS)
    set(BOOST_TEST_SOURCES ${BOOST_TEST_DEFAULT_ARGS})
  else (BOOST_TEST_DEFAULT_ARGS)
    set(BOOST_TEST_SOURCES "${testname}.cpp")
  endif (BOOST_TEST_DEFAULT_ARGS)

  #message("Sources: ${BOOST_TEST_SOURCES}")
  if (BOOST_TEST_RUN)
    set(BOOST_TEST_TAG "run")
  elseif(BOOST_TEST_COMPILE)
    set(BOOST_TEST_TAG "compile")
  elseif(BOOST_TEST_LINK)
    set(BOOST_TEST_TAG "link")
  endif(BOOST_TEST_RUN)

  if (BOOST_TEST_FAIL)
    set(BOOST_TEST_TAG ${BOOST_TEST_TAG}-fail)
  endif(BOOST_TEST_FAIL)

  set(BOOST_TEST_TESTNAME "${PROJECT_NAME}-${testname}-${BOOST_TEST_TAG}")
  #message("testname: ${BOOST_TEST_TESTNAME}")
  # If testing is turned off, this test is not okay
  if (NOT BUILD_TESTING)
    set(BOOST_TEST_OKAY FALSE)
  endif(NOT BUILD_TESTING)

endmacro(boost_test_parse_args)

# This macro creates a Boost regression test that will be executed. If
# the test can be built, executed, and exits with a return code of
# zero, it will be considered to have passed.
#
#   boost_test_run(testname 
#                  [source1 source2 ...]
#                  [ARGS arg1 arg2... ]
#                  [COMPILE_FLAGS compileflags]
#                  [LINK_FLAGS linkflags]
#                  [LINK_LIBS linklibs]
#                  [DEPENDS libdepend1 libdepend2 ...]
#                  [EXTRA_OPTIONS option1 option2 ...])
#
# testname is the name of the test. source1, source2, etc. are the
# source files that will be built and linked into the test
# executable. If no source files are provided, the file "testname.cpp"
# will be used instead.
#
# There are several optional arguments to control how the regression
# test is built and executed:
#
#   ARGS: Provides additional arguments that will be passed to the
#   test executable when it is run.
#
#   COMPILE_FLAGS: Provides additional compilation flags that will be
#   used when building this test. For example, one might want to add
#   "-DBOOST_SIGNALS_ASSERT=1" to turn on assertions within the library.
#
#   LINK_FLAGS: Provides additional flags that will be passed to the
#   linker when linking the test excecutable. This option should not
#   be used to link in additional libraries; see LINK_LIBS and
#   DEPENDS.
#
#   LINK_LIBS: Provides additional libraries against which the test
#   executable will be linked. For example, one might provide "expat"
#   as options to LINK_LIBS, to state that this executable should be
#   linked against the external "expat" library. Use LINK_LIBS for
#   libraries external to Boost; for Boost libraries, use DEPENDS.
#
#   DEPENDS: States that this test executable depends on and links
#   against another Boost library. The argument to DEPENDS should be
#   the name of a particular variant of a Boost library, e.g.,
#   boost_signals-static.
#
#   EXTRA_OPTIONS: Provide extra options that will be passed on to 
#   boost_add_executable.
#
# Example:
#   boost_test_run(signal_test DEPENDS boost_signals)
macro(boost_test_run testname)
  boost_test_parse_args(${testname} ${ARGN} RUN)
  if (BOOST_TEST_OKAY)
#    add_executable(${testname} ${BOOST_TEST_SOURCES})
    boost_add_executable(${testname} ${BOOST_TEST_SOURCES}
      OUTPUT_NAME tests/${PROJECT_NAME}/${testname}
      DEPENDS "${BOOST_TEST_DEPENDS}"
      LINK_LIBS ${BOOST_TEST_LINK_LIBS}
      LINK_FLAGS ${BOOST_TEST_LINK_FLAGS}
      COMPILE_FLAGS ${BOOST_TEST_COMPILE_FLAGS}
      NO_INSTALL 
      ${BOOST_TEST_EXTRA_OPTIONS})
    set (THIS_EXE_OKAY TRUE)
    if (THIS_EXE_OKAY)
      # This target builds and runs the test
      add_custom_target(${BOOST_TEST_TESTNAME})

      file( TO_NATIVE_PATH "${BOOST_TEST_DRIVER}" NATIVE_BOOST_TEST_DRIVER )

      set(THIS_TEST_PREFIX_ARGS
          ${PYTHON_EXECUTABLE} ${NATIVE_BOOST_TEST_DRIVER} 
          ${CMAKE_CURRENT_BINARY_DIR} ${BOOST_TEST_TAG} ${testname} 
	       )

      add_dependencies(${BOOST_TEST_TESTNAME}
    	${testname}
    	)

      add_custom_command(TARGET ${BOOST_TEST_TESTNAME}
                        POST_BUILD
                        COMMAND 
                        ${THIS_TEST_PREFIX_ARGS} 
                        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_CFG_INTDIR}/tests/${PROJECT_NAME}/${testname}
                        ${BOOST_TEST_ARGS}
                        COMMENT "Running ${testname} in project ${PROJECT_NAME}"
                        )

      add_dependencies(${BOOST_TEST_TESTNAME}
    	${PROJECT_NAME}-${testname}
    	)

      add_dependencies(${PROJECT_NAME}-test
    	${BOOST_TEST_TESTNAME}
    	)

    endif(THIS_EXE_OKAY)
  endif (BOOST_TEST_OKAY)
endmacro(boost_test_run)

# 
# This macro creates a boost regression test that will be run but is
# expected to fail (exit with nonzero return code).
# See boost_test_run()
# 
macro(boost_test_run_fail testname)
  boost_test_run(${testname} ${ARGN} FAIL)
endmacro(boost_test_run_fail)


# This macro creates a Boost regression test that will be compiled,
# but not linked or executed. If the test can be compiled with no
# failures, the test passes.
#
#   boost_test_compile(testname 
#                      [source1]
#                      [COMPILE_FLAGS compileflags])
#
# testname is the name of the test. source1 is the name of the source
# file that will be built. If no source file is provided, the file
# "testname.cpp" will be used instead.
#
# The COMPILE_FLAGS argument provides additional arguments that will
# be passed to the compiler when building this test.

# Example:
#   boost_test_compile(advance)
macro(boost_test_compile testname)
  boost_test_parse_args(${testname} ${ARGN} COMPILE)


  if (BOOST_TEST_OKAY)
    # Determine the include directories to pass along to the underlying
    # project.
    # works but not great
    get_directory_property(BOOST_TEST_INCLUDE_DIRS INCLUDE_DIRECTORIES)
    set(BOOST_TEST_INCLUDES "")
    foreach(DIR ${BOOST_TEST_INCLUDE_DIRS})
      set(BOOST_TEST_INCLUDES "${BOOST_TEST_INCLUDES};-I${DIR}")
    endforeach(DIR ${BOOST_TEST_INCLUDE_DIRS})

    set(THIS_TEST_PREFIX_ARGS
      ${PYTHON_EXECUTABLE} ${BOOST_TEST_DRIVER} ${CMAKE_CURRENT_BINARY_DIR} ${BOOST_TEST_TAG} ${testname} 
      )
  
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_TEST_TESTNAME}.${CMAKE_CXX_OUTPUT_EXTENSION}
      COMMAND
      ${THIS_TEST_PREFIX_ARGS}
      ${CMAKE_CXX_COMPILER} 
      ${BOOST_TEST_COMPILE_FLAGS} 
      ${BOOST_TEST_INCLUDES}
      -c ${BOOST_TEST_SOURCES}
      -o ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_TEST_TESTNAME}${CMAKE_CXX_OUTPUT_EXTENSION}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS ${BOOST_TEST_SOURCES}
      COMMENT "Running ${testname} in project ${PROJECT_NAME}"
      )

    add_custom_target(${BOOST_TEST_TESTNAME}
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${BOOST_TEST_TESTNAME}.${CMAKE_CXX_OUTPUT_EXTENSION}
      )

    add_dependencies(${PROJECT_NAME}-test
      ${BOOST_TEST_TESTNAME}
      )

  endif(BOOST_TEST_OKAY)
endmacro(boost_test_compile)

#
# This macro creates a Boost regression test that is expected to 
# *fail* to compile.   See boost_test_compile()
#
macro(boost_test_compile_fail testname)
  boost_test_compile(${testname} ${ARGN} FAIL)
endmacro(boost_test_compile_fail)




#
# boost_test_link:
#
# Under construction.
#
macro(boost_test_link testname)
  boost_test_parse_args(${testname} ${ARGN} LINK)
  if(BOOST_TEST_OKAY)

    set(THIS_TEST_PREFIX_ARGS
#      ${PYTHON_EXECUTABLE} ${BOOST_TEST_DRIVER} ${CMAKE_CURRENT_BINARY_DIR} test_link ${testname} 
      )
    
    #
    #  FIXME:  no ctest.
    #
    add_custom_target(TARGET ${BOOST_TEST_TESTNAME}
      COMMAND /link/tests/are/failing/at/the/moment
      COMMENT "Link test ${testname} in ${PROJECT_NAME} is failing."
      )

    # POST_BUILD
    # COMMAND 
    # ${THIS_TEST_PREFIX_ARGS}
    # ${CMAKE_CTEST_COMMAND}
    # --build-and-test
    # ${Boost_SOURCE_DIR}/tools/build/CMake/LinkTest
    # ${Boost_BINARY_DIR}/tools/build/CMake/LinkTest
    # --build-generator \\"${CMAKE_GENERATOR}\\"
    # --build-makeprogram \\"${MAKEPROGRAM}\\"
    # --build-project LinkTest
    # --build-options -DSOURCE=${CMAKE_CURRENT_SOURCE_DIR}/${BOOST_TEST_SOURCES} -DINCLUDES=${Boost_SOURCE_DIR} -DCOMPILE_FLAGS=\\"${BOOST_TEST_COMPILE_FLAGS}\\"
    # COMMENT "Running ${testname} (link) in project ${PROJECT_NAME}"
    # )

    add_dependencies(${PROJECT_NAME}-test
      ${BOOST_TEST_TESTNAME}
      )
    
  endif(BOOST_TEST_OKAY)
endmacro(boost_test_link)

