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
      filesystem serialization system program_options regex
      thread date_time chrono timer iostreams unit_test_framework)

  if(ALPS_ENABLE_MPI)
    list(APPEND _alps_boost_components mpi)
  endif()

  # Find the non-Python components first (all required).
  find_package(Boost REQUIRED COMPONENTS ${_alps_boost_components})

  # Python component: library naming varies by Boost/distro version.
  # Try python<major><minor> (e.g. python311), then python3, then python.
  if(ALPS_HAVE_PYTHON)
    set(_alps_python_component "")
    foreach(_pycomp "python${PYVER}" "python3" "python")
      find_package(Boost QUIET COMPONENTS ${_pycomp})
      if(Boost_${_pycomp}_FOUND)
        set(_alps_python_component ${_pycomp})
        break()
      endif()
    endforeach()

    if(_alps_python_component)
      message(STATUS "Found system Boost.Python component: ${_alps_python_component}")
      # Merge the python library into Boost_LIBRARIES if not already there.
      list(APPEND Boost_LIBRARIES Boost::${_alps_python_component})
    else()
      message(WARNING
        "System Boost.Python library not found (tried python${PYVER}, python3, python). "
        "Python bindings will be disabled.")
      set(ALPS_HAVE_PYTHON OFF)
      set(BUILD_BOOST_PYTHON OFF)
    endif()
  endif()

  # Set ALPS_HAVE_BOOST_NUMPY for Boost >= 1.63 (when boost::python::numpy
  # was introduced).
  if(Boost_VERSION_STRING VERSION_GREATER_EQUAL "1.63.0")
    set(ALPS_HAVE_BOOST_NUMPY ON)
  endif()

  # Scenario 3: system Boost 1.63-1.86 + NumPy >= 2.0.
  # boost::python::numpy in these versions uses deprecated NumPy C API
  # removed in NumPy 2.0.  Fall back to boost::python::numeric::array,
  # which uses only the stable NumPy C API.
  if(ALPS_HAVE_BOOST_NUMPY AND ALPS_HAVE_PYTHON)
    EXEC_PYTHON_SCRIPT("import numpy; print(numpy.__version__)" _alps_numpy_ver)
    message(STATUS "NumPy version: ${_alps_numpy_ver}")
    if(_alps_numpy_ver VERSION_GREATER_EQUAL "2.0.0" AND
       Boost_VERSION_STRING VERSION_LESS "1.87.0")
      message(WARNING
        "System Boost ${Boost_VERSION_STRING} does not support NumPy >= 2.0 "
        "(requires Boost >= 1.87). "
        "Falling back to boost::python::numeric::array. "
        "Upgrade system Boost to >= 1.87 to silence this warning.")
      set(ALPS_HAVE_BOOST_NUMPY OFF)
    endif()
  endif()

  # Align Boost_INCLUDE_DIR (singular) used elsewhere in the build.
  set(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIRS})

  message(STATUS "Using system Boost ${Boost_VERSION_STRING}")
  message(STATUS "  includes:           ${Boost_INCLUDE_DIRS}")
  message(STATUS "  libraries:          ${Boost_LIBRARIES}")
  message(STATUS "  ALPS_HAVE_BOOST_NUMPY: ${ALPS_HAVE_BOOST_NUMPY}")

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
