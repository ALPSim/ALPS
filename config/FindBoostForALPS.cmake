# - FindBoostForALPS.cmake
# Find Boost precompiled libraries or Boost source tree for ALPS
#

#  Copyright Ryo IGARASHI 2010, 2011, 2013.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

# ALPS_USE_SYSTEM_BOOST: when ON, use a system-installed precompiled Boost
# instead of building from source.  The source build remains the default.
option(ALPS_USE_SYSTEM_BOOST "Use system-installed Boost instead of building from source" OFF)

if(ALPS_USE_SYSTEM_BOOST)

  # -----------------------------------------------------------------------
  # System Boost path
  # -----------------------------------------------------------------------

  set(_alps_boost_components
      filesystem serialization program_options regex
      thread date_time chrono timer iostreams unit_test_framework)

  if(ALPS_ENABLE_MPI)
    list(APPEND _alps_boost_components mpi)
  endif()

  # Force CMake module mode (FindBoost.cmake) so component names like
  # unit_test_framework work regardless of whether the system Boost ships
  # a BoostConfig.cmake (config mode).  Config mode uses different component
  # names and often conflicts with BOOST_ROOT overrides.
  set(Boost_NO_BOOST_CMAKE TRUE)
  # CMP0167 (CMake 3.30+) warns that the FindBoost module is deprecated in
  # favour of Boost's own config file.  We intentionally use module mode here,
  # so acknowledge the policy to suppress the warning.
  if(POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
  endif()

  # Find the non-Python components first (all required).
  # Note: boost::system is intentionally absent — it became header-only in
  # Boost 1.69 (no library file to link).  Boost < 1.69 is not supported
  # in this mode; the version check below enforces that.
  find_package(Boost REQUIRED COMPONENTS ${_alps_boost_components})

  # Enforce the minimum: for Boost < 1.69 boost::system is a compiled
  # library required by filesystem, and omitting it from the component list
  # would produce a mis-linked binary.
  if(Boost_VERSION_STRING VERSION_LESS "1.69.0")
    message(FATAL_ERROR
      "ALPS_USE_SYSTEM_BOOST requires Boost >= 1.69 "
      "(found ${Boost_VERSION_STRING}). "
      "Upgrade the system Boost installation or disable ALPS_USE_SYSTEM_BOOST.")
  endif()

  # Align Boost_INCLUDE_DIR (singular) used elsewhere in the build.
  set(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIRS})

  # Add the Boost lib directory to the linker search path so that @rpath
  # transitive dependencies of Boost dylibs (e.g. libboost_graph,
  # libboost_container) are found at link time on macOS.
  if(Boost_LIBRARY_DIRS)
    link_directories(${Boost_LIBRARY_DIRS})
  endif()

  message(STATUS "Using system Boost ${Boost_VERSION_STRING}")
  message(STATUS "  includes:           ${Boost_INCLUDE_DIRS}")
  message(STATUS "  libraries:          ${Boost_LIBRARIES}")
  return()

endif() # ALPS_USE_SYSTEM_BOOST

# -----------------------------------------------------------------------
# Source-build path (original behaviour)
# -----------------------------------------------------------------------

# Since Boost_ROOT_DIR is used for setting Boost source directory,
# we use precompiled Boost libraries only when Boost_ROOT_DIR is not set.

if (NOT Boost_SRC_DIR)
  message(STATUS "Environment variable Boost_SRC_DIR has not been set.")
  message(STATUS "Downloading the most recent Boost release...")

  include(FetchContent)

  # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
  if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    cmake_policy(SET CMP0135 NEW)
  endif()

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
