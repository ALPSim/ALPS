# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#if (NOT Boost_ROOT_DIR)
#  find_package(Boost 1.42.0 COMPONENTS date_time filesystem program_options system python thread)
#endif(NOT BOOST_ROOT_DIR)
message(STATUS "Looking for Boost Source")
find_package(BoostSrc)
