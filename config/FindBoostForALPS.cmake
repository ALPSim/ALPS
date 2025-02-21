# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010, 2011, 2013.
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

# Since Boost_ROOT_DIR is used for setting Boost source directory,
# we use precompiled Boost libraries only when Boost_ROOT_DIR is not set.

if (NOT Boost_SRC_DIR)
  message(STATUS "Environment variable Boost_SRC_DIR has not been set.")
  message(STATUS "Downloading the most recent Boost release...")

  include(FetchContent)

  FetchContent_Declare(
    boost_src
    URL      https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz
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
