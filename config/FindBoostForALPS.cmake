# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010, 2011, 2013.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Since Boost_ROOT_DIR is used for setting Boost source directory,
# we use precompiled Boost libraries only when Boost_ROOT_DIR is not set.

if (NOT Boost_SRC_DIR)
  message(STATUS "Environment variable Boost_SRC_DIR has not been set.")
  message(STATUS "Downloading the most recent Boost release...")

  include(FetchContent)

  FetchContent_Declare(
    boost_src
    URL      https://boostorg.jfrog.io/artifactory/main/release/1.87.0/source/boost_1_87_0_rc1.tar.gz
    URL_HASH SHA256=f55c340aa49763b1925ccf02b2e83f35fdcf634c9d5164a2acb87540173c741d
    EXCLUDE_FROM_ALL
  )

  FetchContent_MakeAvailable(boost_src)

  message(STATUS "Boost sources are in ${boost_src_SOURCE_DIR}")

#  set(Boost_FOUND TRUE)
  set(Boost_ROOT_DIR ${boost_src_SOURCE_DIR})

else(NOT Boost_SRC_DIR)
  set(Boost_ROOT_DIR ${Boost_SRC_DIR})

endif(NOT Boost_SRC_DIR)

# Boost_FOUND is set only when FindBoost.cmake succeeds.
# if not, build Boost libraries from source.
if (NOT Boost_FOUND)
  message(STATUS "Looking for Boost Source")
  find_package(BoostSrc)
endif(NOT Boost_FOUND)
