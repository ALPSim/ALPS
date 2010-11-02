# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Since Boost_ROOT_DIR is used for setting Boost source directory,
# we use precompiled Boost libraries only when Boost_ROOT_DIR is not set.
if (NOT Boost_ROOT_DIR)
  message(STATUS "Looking for precompiled Boost libraries")
  # Boost 1.42 is not (yet) supported version number on CMake 1.8.0.
  set(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0")
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
  find_package(Boost 1.42.0 COMPONENTS mpi)
  if (Boost_FOUND)
    find_package(Boost 1.42.0 COMPONENTS date_time filesystem program_options python regex system serialization thread mpi)
  else (Boost_FOUND)
    find_package(Boost 1.42.0 COMPONENTS date_time filesystem program_options python regex system serialization thread)
  endif (Boost_FOUND)
  if (APPLE AND NOT Boost_FOUND)
    find_package(boost COMPONENTS date_time-mt filesystem-mt program_options-mt python-mt regex-mt system-mt serialization-mt thread-mt mpi-mt)
  endif(APPLE AND NOT Boost_FOUND)
endif(NOT Boost_ROOT_DIR)

# Boost_FOUND is set only when FindBoost.cmake succeeds.
if (NOT Boost_FOUND)
  message(STATUS "Looking for Boost Source")
  find_package(BoostSrc)
endif(NOT Boost_FOUND)
