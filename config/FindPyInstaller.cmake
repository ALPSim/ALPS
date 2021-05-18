# Find PyInstaller (PI)
# Will find the path to Makespec.py and Build.py

# Look for Python:
find_package(PythonLibs REQUIRED)

if(NOT PyInstaller_PATH)
  find_path(PyInstaller_PATH makespec.py
            $ENV{HOME}/pyinstaller
            $ENV{HOME}/src/pyinstaller
            $ENV{HOME}/ALPS/src/pyinstaller
            $ENV{HOME}/pyinstaller-1.3
            $ENV{HOME}/src/pyinstaller-1.3
            $ENV{HOME}/ALPS/src/pyinstaller-1.3
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/pyinstaller"
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/src/pyinstaller"
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/ALPS/src/pyinstaller"
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/pyinstaller-1.3"
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/src/pyinstaller-1.3"
            "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/My Documents/ALPS/src/pyinstaller-1.3"
            DOC "Path to the pyinstaller directory (where to find makespec.py)")
  if(PyInstaller_PATH)
    set(PyInstaller_PATH ${PyInstaller_PATH} CACHE PATH "Path to PyInstaller")
    set(PyInstaller_FOUND 1 CACHE BOOL "Found PyInstaller")
    message(STATUS "Found PyInstaller: ${PyInstaller_PATH}")
  else(PyInstaller_PATH)
    set(PyInstaller_PATH "" CACHE PATH "Path to PyInstaller")
    set(PyInstaller_FOUND 0 CACHE BOOL "Found PyInstaller")
    message(STATUS "PyInstaller: not found")
  endif(PyInstaller_PATH)
endif(NOT PyInstaller_PATH)

if(PyInstaller_PATH)
  set(PI_MAKESPEC ${PyInstaller_PATH}/makespec.py)
  set(PI_BUILD ${PyInstaller_PATH}/Build.py)
endif(PyInstaller_PATH)

mark_as_advanced(PI_MAKESPEC PI_BUILD)
