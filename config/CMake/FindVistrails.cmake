# look for Vistrails and set the following variables:

# Vistrails_FOUND
# VISTRAILS_APP_NAME
# VISTRAILS_APP_DIR
# VISTRAILS_PACKAGE_DIR
# VISTRAILS_PYTHON_EXTENSION_DIR
# VISTRAILS_LIB_DIR
# VISTRAILS_DYLIB_DIR
# VISTRAILS_PYTHON_INTERPRETER

if(WIN32 AND NOT UNIX)
  set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the Vistrails application")
  set(VISTRAILS_PACKAGE_DIR "packages")
  set(VISTRAILS_PYTHON_EXTENSION_DIR Python25/Lib/site-packages )
  set(VISTRAILS_LIB_DIR "lib")
  set(VISTRAILS_PYTHON_EXTENSION_DIR "Python25/DLLs")
  set(VISTRAILS_DYLIB_DIR "")
else(WIN32 AND NOT UNIX)
  if(APPLE)	
    set(VISTRAILS_APP_NAME "Vistrails.app" CACHE STRING "Name of the Vistrails application")
    find_path(VISTRAILS_APP_DIR ${VISTRAILS_APP_NAME} "/Applications/Vistrails/")
    if(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "TRUE")
    else(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "FALSE")
    endif(VISTRAILS_APP_DIR)
    set(VISTRAILS_PACKAGE_DIR "Contents/Resources/lib/python2.5/packages")
    set(VISTRAILS_PYTHON_EXTENSION_DIR "Contents/Resources/lib/python2.5")
    set(VISTRAILS_LIB_DIR "Contents/Resources/lib")
    set(VISTRAILS_DYLIB_DIR "Contents/Frameworks")
    set(VISTRAILS_PYTHON_INTERPRETER "Contents/MacOS/python")
  else(APPLE)	
    set(VISTRAILS_APP_NAME "vistrails" CACHE STRING "Name of the Vistrails application")
    find_path(VISTRAILS_APP_DIR ${VISTRAILS_APP_NAME} ${PATH})
    if(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "TRUE")
    else(VISTRAILS_APP_DIR)
      set(VISTRAILS_FOUND "FALSE")
    endif(VISTRAILS_APP_DIR)
    set(VISTRAILS_PACKAGE_DIR "vistrails/packages")
    set(VISTRAILS_PYTHON_EXTENSION_DIR "/usr/lib/python2.5/site-packages")
    set(VISTRAILS_LIB_DIR "vistrails/lib")
    set(VISTRAILS_DYLIB_DIR "lib")
  endif(APPLE)
endif(WIN32 AND NOT UNIX)
