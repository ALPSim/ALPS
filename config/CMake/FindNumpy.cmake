# Find PyInstaller (PI)
# Will find the path to Makespec.py and Build.py

# Look for Python:
find_package(PythonLibs REQUIRED)

if (PYTHONLIBS_FOUND AND NOT PYTHON_NUMPY_INCLUDE_DIR)
    set(SEARCH_DIRS ${PYTHON_INCLUDE_DIR}/../Versions/Current/Extras/lib/python/numpy/core/include)
    if(VISTRAILS_APP_NAME)
        set(SEARCH_DIRS ${SEARCH_DIRS} 
           ${CMAKE_INSTALL_PREFIX}/${VISTRAILS_APP_NAME}/${VISTRAILS_PYTHON_EXTENSION_DIR}/numpy/core/include
        )
    endif(VISTRAILS_APP_NAME)
    find_path(NUMPY_INCLUDE_DIR numpy/ndarrayobject.h ${SEARCH_DIRS})
            
    if(NUMPY_INCLUDE_DIR)
        set(PYTHON_NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE PATH "Path to NUMPY C API")
    endif(NUMPY_INCLUDE_DIR)

endif (PYTHONLIBS_FOUND AND NOT PYTHON_NUMPY_INCLUDE_DIR)

mark_as_advanced(PYTHON_NUMPY_INCLUDE_DIR)
