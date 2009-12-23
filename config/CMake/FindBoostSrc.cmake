# - FindBoostSrc
# Find Boost source tree
#
#
#

if(BOOST_ROOT)
  find_path(Boost_ROOT_DIR "libs/filesystem/src/path.cpp" ${BOOST_ROOT})
else(BOOST_ROOT)
  set(DIR0
    $ENV{HOME} $ENV{HOME}/src $ENV{HOME}/ALPS/src /usr/local /usr/local/src
    "$ENV{HOMEDRIVE}/Program Files"
    "$ENV{HOMEDRIVE}$ENV{HOMEPATH}"
    "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/src"
    "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/ALPS/src")
  set(DIR1 boost boostsrc boost_1_41_0 boost_1_40_0 boost_1_39_0 boost_1_38_0 boostsrc_1_38_0)
  set(_boost_SEARCH_PATH "")
  foreach(D0 ${DIR0})
    foreach(D1 ${DIR1})
      set(_boost_SEARCH_PATH ${_boost_SEARCH_PATH} ${D0}/${D1})
    endforeach(D1)
  endforeach(D0)
  find_path(Boost_ROOT_DIR "libs/filesystem/src/path.cpp" ${_boost_SEARCH_PATH})
endif(BOOST_ROOT)

if(Boost_ROOT_DIR)
  set(Boost_INCLUDE_DIR ${Boost_ROOT_DIR} CACHE PATH "Boost Include Directory")
endif(Boost_ROOT_DIR)

# check Boost version (from Modules/FindBoost.cmake)
if(Boost_INCLUDE_DIR)
  set(_boost_VERSION 0)
  set(_boost_LIB_VERSION "")
  file(READ "${Boost_INCLUDE_DIR}/boost/version.hpp" _boost_VERSION_HPP_CONTENTS)
  string(REGEX REPLACE ".*#define BOOST_VERSION ([0-9]+).*" "\\1" _boost_VERSION "${_boost_VERSION_HPP_CONTENTS}")
  string(REGEX REPLACE ".*#define BOOST_LIB_VERSION \"([0-9_]+)\".*" "\\1" _boost_LIB_VERSION "${_boost_VERSION_HPP_CONTENTS}")
  set(Boost_LIB_VERSION ${_boost_LIB_VERSION} CACHE INTERNAL "The library version string for boost libraries")
  set(Boost_VERSION ${_boost_VERSION} CACHE INTERNAL "The version number for boost libraries")
  if(NOT "${Boost_VERSION}" STREQUAL "0")
    MATH(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
    MATH(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
    MATH(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
  endif(NOT "${Boost_VERSION}" STREQUAL "0")
endif(Boost_INCLUDE_DIR)

if(Boost_ROOT_DIR)
  message(STATUS "Found Boost Source: ${Boost_ROOT_DIR}")
else(Boost_ROOT_DIR)
  message(FATAL_ERROR "Boost Source not Found")
endif(Boost_ROOT_DIR)

set(BUILD_BOOST_DATE_TIME TRUE)
set(BUILD_BOOST_FILESYSTEM TRUE)
set(BUILD_BOOST_PROGRAM_OPTIONS TRUE)
set(BUILD_BOOST_REGEX TRUE)
set(BUILD_BOOST_SERIALIZATION TRUE)
set(BUILD_BOOST_SYSTEM TRUE)
set(BUILD_BOOST_PYTHON TRUE)

# Avoid auto link of Boost library
add_definitions(-DBOOST_ALL_NO_LIB=1)
if(BUILD_SHARED_LIBS)
  add_definitions(-DBOOST_ALL_DYN_LINK=1)
else(BUILD_SHARED_LIBS)
  add_definitions(-DBOOST_ALL_STATIC_LINK=1)
endif(BUILD_SHARED_LIBS)

mark_as_advanced(Boost_INCLUDE_DIR)
