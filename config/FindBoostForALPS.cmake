# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010, 2011, 2013.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Since Boost_ROOT_DIR is used for setting Boost source directory,
# we use precompiled Boost libraries only when Boost_ROOT_DIR is not set.
if (NOT Boost_ROOT_DIR)

  # Old version of CMake (< 2.8.8) do not support for new Boost version
  set(Boost_ADDITIONAL_VERSIONS "1.56.0" "1.56" "1.55.0" "1.55" "1.54.0" "1.54" "1.53.0" "1.53" "1.52.0" "1.52")
  # Debian and Ubuntu packages are multithreaded by default but not on MacOSX.
  # Explicitly set multithread flag.
  set(Boost_USE_MULTITHREADED ON)
  # If you want to use static libraries, uncomment this option.
  #  set(Boost_USE_STATIC_LIBS ON)
  # Debug flag for FindBoost.cmake
  #  set(Boost_DEBUG TRUE)

  message(STATUS "Looking for precompiled Boost libraries (version >= 1.49)")

  # We do not require Boost.MPI, therefore check whether Boost.MPI exists
  # before actual find_package for Boost.
  # - Ubuntu 10.10 does not have Boost.MPI package.
  set(MPI)
  find_package(Boost 1.49.0 COMPONENTS mpi)
  if (Boost_FOUND)
    set(MPI mpi)
  endif (Boost_FOUND)
  # Linking to Boost.test Unit Test Framework is optional
  set(UTF)
  if (ALPS_LINK_BOOST_TEST)
    find_package(Boost 1.49.0 COMPONENTS unit_test_framework)
    if (Boost_FOUND)
      set(UTF unit_test_framework)
      # Unset ALPS_INSTALL_BOOST_TEST
      # since Boost.Test dynamic library is already installed
      unset(ALPS_INSTALL_BOOST_TEST)
      # Add definitions needed to link Boost.test Unit Test Framework.
      add_definitions(-DBOOST_TEST_DYN_LINK)
      add_definitions(-DALPS_LINK_BOOST_TEST)
    endif (Boost_FOUND)
  endif(ALPS_LINK_BOOST_TEST)
  # Linking to Boost.timer (v2) is not used (and thus optional) as of ALPS 2.1.
  # However, Boost.timer (v2) is built and linked when boost is built together.
  # Therefore, When Boost.timer (v2) is found, link ALPS against it.
  set(TIMER)
  find_package(Boost 1.48.0 COMPONENTS timer)
  if (Boost_FOUND)
    set(TIMER timer)
  endif (Boost_FOUND)
  set(NUMPY)
  set(PYTHON python)
  set(ALPS_HAVE_BOOST_NUMPY OFF)
  find_package(Boost 1.48.0 COMPONENTS numpy)
  if (Boost_FOUND)
    set(NUMPY numpy)
    set(PYTHON python)
    set(ALPS_HAVE_BOOST_NUMPY ON)
  endif (Boost_FOUND)
  find_package(Boost 1.48.0 COMPONENTS numpy27)
  if (Boost_FOUND)
    set(NUMPY numpy27)
    set(PYTHON python27)
    set(ALPS_HAVE_BOOST_NUMPY ON)
  endif (Boost_FOUND)
  # The final library finding for precompiled Boost.
  if(ALPS_BUILD_PYTHON)
    find_package(Boost 1.49.0 COMPONENTS chrono date_time filesystem iostreams program_options ${PYTHON} regex system serialization thread ${MPI} ${UTF} ${TIMER} ${NUMPY})
  else(ALPS_BUILD_PYTHON)
    find_package(Boost 1.49.0 COMPONENTS chrono date_time filesystem iostreams program_options regex system serialization thread ${MPI} ${UTF} ${TIMER})
  endif(ALPS_BUILD_PYTHON)
  # Unset local variable
  unset(MPI)
  unset(UTF)
  unset(TIMER)
  unset(NUMPY)
  unset(PYTHON)

endif(NOT Boost_ROOT_DIR)

# Boost_FOUND is set only when FindBoost.cmake succeeds.
# if not, build Boost libraries from source.
if (NOT Boost_FOUND)
  message(STATUS "Looking for Boost Source")
  find_package(BoostSrc)
endif(NOT Boost_FOUND)
