#  Copyright ALPS collaboration 2026.
#   Distributed under the MIT licence; see LICENSE.txt.
#
# Single source of truth for the ALPS release version.
#
# ALPS_VERSION.txt holds the numeric release core (MAJOR.MINOR.PATCH) and
# nothing else. Two constraints force that:
#
#   * project(VERSION ...) rejects anything non-numeric, so "2.4.0-beta.1"
#     fails to configure.
#   * find_package() version matching and the library SOVERSION have no notion
#     of prerelease ordering.
#
# A prerelease label such as "beta.2" therefore lives in ALPS_VERSION_PRERELEASE
# (see the version block in the top-level CMakeLists.txt), where it affects the
# display string only -- never the numeric version used for ABI and
# find_package() decisions.
#
# This file is included by full path before project(), so it cannot rely on
# CMAKE_MODULE_PATH or PROJECT_SOURCE_DIR.

set(_alps_version_file "${CMAKE_CURRENT_LIST_DIR}/../ALPS_VERSION.txt")

if(NOT EXISTS "${_alps_version_file}")
  message(FATAL_ERROR "Cannot read the ALPS version file: ${_alps_version_file}")
endif()

file(STRINGS "${_alps_version_file}" ALPS_VERSION_CORE LIMIT_COUNT 1)
string(STRIP "${ALPS_VERSION_CORE}" ALPS_VERSION_CORE)

# Fail loudly here rather than letting project() emit "VERSION format invalid",
# which gives no hint about which file is at fault.
if(NOT ALPS_VERSION_CORE MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
  message(FATAL_ERROR
    "${_alps_version_file} must contain exactly MAJOR.MINOR.PATCH, but reads "
    "'${ALPS_VERSION_CORE}'. Prerelease labels belong in "
    "ALPS_VERSION_PRERELEASE, and the leading 'v' of a release tag is not "
    "part of the version.")
endif()

unset(_alps_version_file)
