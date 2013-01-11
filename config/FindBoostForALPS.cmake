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
  set(Boost_ADDITIONAL_VERSIONS "1.56.0" "1.56" "1.55.0" "1.55" "1.54.0" "1.54" "1.53.0" "1.53" "1.52.0" "1.52" "1.51.0" "1.51" "1.50.0" "1.50" "1.49.0" "1.49" "1.48.0" "1.48" "1.47.0" "1.47")
  # Debian and Ubuntu packages are multithreaded.
  # Commented out because default is ON.
  #  set(Boost_USE_MULTITHREADED ON)
  # If you want to use static libraries, uncomment this option.
  #  set(Boost_USE_STATIC_LIBS ON)
  # Debug flag for FindBoost.cmake
  #  set(Boost_DEBUG TRUE)

  # We do not require Boost.MPI, therefore check whether Boost.MPI exists
  # before actual find_package for Boost.
  # - Ubuntu 10.10 does not have Boost.MPI package.
    # Boost 1.47.0 or later is not yet supported by FindBoost.cmake (CMake 2.8.5).
    message(STATUS "Looking for precompiled Boost libraries (version >= 1.48)")

    find_package(Boost 1.48.0 COMPONENTS mpi)
    if (ALPS_LINK_BOOST_TEST)
      set(UTF "unit_test_framework")
    else (ALPS_LINK_BOOST_TEST)
      set(UTF)
    endif(ALPS_LINK_BOOST_TEST)
    if (Boost_FOUND)
      find_package(Boost 1.48.0 COMPONENTS chrono date_time filesystem program_options python regex system serialization thread mpi ${UTF})
    else (Boost_FOUND)
      find_package(Boost 1.48.0 COMPONENTS chrono date_time filesystem program_options python regex system serialization thread ${UTF})
    endif (Boost_FOUND)
    if (APPLE AND NOT Boost_FOUND)
      find_package(boost COMPONENTS chrono-mt date_time-mt filesystem-mt program_options-mt python-mt regex-mt system-mt serialization-mt thread-mt mpi-mt ${UTF}-mt)
    endif(APPLE AND NOT Boost_FOUND)
    if (UTF AND Boost_FOUND)
      # Unset ALPS_INSTALL_BOOST_TEST
      # since Boost.Test dynamic library is already installed
      unset(ALPS_INSTALL_BOOST_TEST)
      # Add definitions needed to link Boost.test Unit Test Framework.
      add_definitions(-DBOOST_TEST_DYN_LINK)
      add_definitions(-DALPS_LINK_BOOST_TEST)
    endif (UTF AND Boost_FOUND)
    # Unset local variable
    unset(UTF)

endif(NOT Boost_ROOT_DIR)

# Boost_FOUND is set only when FindBoost.cmake succeeds.
if (NOT Boost_FOUND)
  message(STATUS "Looking for Boost Source")
  find_package(BoostSrc)
endif(NOT Boost_FOUND)
