# Find PyInstaller (PI)
# Will find the path to Makespec.py and Build.py

# Look for Python:

EXECUTE_PROCESS(COMMAND ${PARCOLLET_INTERPRETER} -c "import distutils ; from distutils.sysconfig import * ; print distutils.sysconfig.get_python_inc()" OUTPUT_VARIABLE PARCOLLET_INCLUDE_DIR)
SET(PARCOLLET_LIBRARIES -Lhere -lthat)
MESSAGE(STATUS "Found it in: ${PARCOLLET_INCLUDE_DIR}")
mark_as_advanced(PARCOLLET_INCLUDE_DIR)
