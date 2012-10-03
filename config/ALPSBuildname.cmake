#  Copyright Michele Dolfi 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# Setting additional information in the CTest build name:
# - NGS enabled?
# - System Boost or path to sources?
# - Version of boost


# Add Boost settings
if(Boost_ROOT_DIR)
  set(BUILD_NAME_SYSTEM_NAME "${CMAKE_SYSTEM_NAME} - Boost Src")
else()
  set(BUILD_NAME_SYSTEM_NAME "${CMAKE_SYSTEM_NAME} - Boost Built")
endif(Boost_ROOT_DIR)

set(BUILD_NAME_SYSTEM_NAME "${BUILD_NAME_SYSTEM_NAME} v${Boost_MAJOR_VERSION}_${Boost_MINOR_VERSION}_${Boost_SUBMINOR_VERSION}")
