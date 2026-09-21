#  Copyright ALPS collaboration 2026.
#   Distributed under the MIT licence; see LICENSE.txt.
# project(VERSION) and package matching need a numeric release core. Keep
# prerelease labels in ALPS_VERSION_PRERELEASE, outside the shared version file.
# This runs before project(), including when ALPS is embedded.
set(_alps_version_file "${CMAKE_CURRENT_LIST_DIR}/../ALPS_VERSION.txt")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_alps_version_file}")
file(READ "${_alps_version_file}" ALPS_VERSION_CORE)
string(STRIP "${ALPS_VERSION_CORE}" ALPS_VERSION_CORE)
if(NOT ALPS_VERSION_CORE MATCHES "^(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)\\.(0|[1-9][0-9]*)$")
  message(FATAL_ERROR
    "${_alps_version_file} must contain exactly MAJOR.MINOR.PATCH, but reads "
    "'${ALPS_VERSION_CORE}'. Prerelease labels belong in "
    "ALPS_VERSION_PRERELEASE, and the leading 'v' of a release tag is not "
    "part of the version.")
endif()
unset(_alps_version_file)
