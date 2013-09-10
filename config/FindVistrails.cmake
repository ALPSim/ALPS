#  Copyright Bela Bauer and Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# look for VisTrails and set the following variables:

# Vistrails_FOUND
# VISTRAILS_APP_NAME
# VISTRAILS_APP_DIR
# VISTRAILS_PACKAGE_DIR
# VISTRAILS_PYTHON_EXTENSION_DIR
# VISTRAILS_LIB_DIR
# VISTRAILS_DYLIB_DIR
# VISTRAILS_PYTHON_INTERPRETER

if(WIN32 AND NOT UNIX)
  set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the VisTrails application")
  if (NOT ALPS_VISTRAILS_32BIT)
    find_path(VISTRAILS_APP_DIR runvistrails.py "$ENV{HOMEDRIVE}/Program Files/VisTrails" "$ENV{HOMEDRIVE}/Program Files/VisTrails" CACHE STRING "Path to the VisTrails directory")
    if (VISTRAILS_APP_DIR)
      set(VISTRAILS_PYTHON_DIR "Python27_64" CACHE STRING "Path to the VisTrails Python directory")
    endif (VISTRAILS_APP_DIR)
  endif (NOT ALPS_VISTRAILS_32BIT)
  if (NOT VISTRAILS_APP_DIR)
    find_path(VISTRAILS_APP_DIR runvistrails.py "$ENV{HOMEDRIVE}/Program Files (x86)/VisTrails" "$ENV{HOMEDRIVE}/Program Files/VisTrails" CACHE STRING "Path to the VisTrails directory")
    if (VISTRAILS_APP_DIR)
      set(VISTRAILS_PYTHON_DIR "Python27" CACHE STRING "Path to the VisTrails Python directory")
    endif (VISTRAILS_APP_DIR)
  endif (NOT VISTRAILS_APP_DIR)
  if(VISTRAILS_APP_DIR)
    set(VISTRAILS_FOUND "TRUE")
  else(VISTRAILS_APP_DIR)
    set(VISTRAILS_FOUND "FALSE")
  endif(VISTRAILS_APP_DIR)
  set(VISTRAILS_PACKAGE_DIR "${VISTRAILS_APP_NAME}/packages")
  set(VISTRAILS_PYTHON_EXTENSION_DIR "${VISTRAILS_PYTHON_DIR}/Lib/site-packages" )
  set(VISTRAILS_LIB_DIR "lib")
  set(VISTRAILS_PYTHONPATH_DIR "")
  set(VISTRAILS_DYLIB_DIR ${VISTRAILS_PYTHON_EXTENSION_DIR}/pyalps)
  set(VISTRAILS_ALTERNATE_APP_DIR "$ENV{HOMEDRIVE}/Program Files (x86)/VisTrails")
  set(VISTRAILS_PYTHON_INTERPRETER "${VISTRAILS_PYTHON_DIR}/python.exe")
else(WIN32 AND NOT UNIX)
  if(APPLE)	
    set(VISTRAILS_APP_NAME "VisTrails.app" CACHE STRING "Name of the VisTrails application")
    find_path(VISTRAILS_APP_DIR ${VISTRAILS_APP_NAME} "/Applications/VisTrails/" CACHE STRING "Path to the VisTrails directory")
    if(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "TRUE")
    else(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "FALSE")
    endif(VISTRAILS_APP_DIR)
      set(VISTRAILS_PACKAGE_DIR "${VISTRAILS_APP_NAME}/Contents/Resources/lib/python2.7/vistrails/packages")
      set(VISTRAILS_PYTHON_EXTENSION_DIR "${VISTRAILS_APP_NAME}/Contents/Resources/lib/python2.7")
    set(VISTRAILS_LIB_DIR "Contents/Resources/lib")
    set(VISTRAILS_PYTHONPATH_DIR "Contents/Resources")
    set(VISTRAILS_DYLIB_DIR "${VISTRAILS_APP_NAME}/Contents/Frameworks")
    set(VISTRAILS_PYTHON_INTERPRETER "Contents/MacOS/python")
  else(APPLE)	
    find_package(Python)
    set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the VisTrails application")
    find_path(VISTRAILS_APP_DIR ${VISTRAILS_APP_NAME} ${PATH})
    if(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "TRUE")
    else(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "FALSE")
    endif(VISTRAILS_APP_DIR)
    set(VISTRAILS_PACKAGE_DIR "${VISTRAILS_APP_NAME}/packages")
    set(VISTRAILS_PYTHON_EXTENSION_DIR ${PYTHON_SITE_PKG} CACHE PATH "Python extensions")
    set(VISTRAILS_LIB_DIR "vistrails/lib")
    set(VISTRAILS_DYLIB_DIR "lib")
  endif(APPLE)
endif(WIN32 AND NOT UNIX)

MARK_AS_ADVANCED( VISTRAILS_APP_NAME )
MARK_AS_ADVANCED( VISTRAILS_PACKAGE_DIR )
MARK_AS_ADVANCED( VISTRAILS_PYTHON_EXTENSION_DIR )
MARK_AS_ADVANCED( VISTRAILS_LIB_DIR )
MARK_AS_ADVANCED( VISTRAILS_DYLIB_DIR )
MARK_AS_ADVANCED( VISTRAILS_PYTHON_INTERPRETER )
MARK_AS_ADVANCED( VISTRAILS_PYTHONPATH_DIR )

