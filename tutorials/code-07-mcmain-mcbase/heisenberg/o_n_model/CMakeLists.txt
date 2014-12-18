cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(ndim_spin NONE)

# find ALPS Library
find_package(ALPS REQUIRED PATHS ${ALPS_ROOT_DIR} $ENV{ALPS_HOME} NO_SYSTEM_ENVIRONMENT_PATH)
message(STATUS "Found ALPS: ${ALPS_ROOT_DIR} (revision: ${ALPS_VERSION})")
include(${ALPS_USE_FILE})

# enable C and C++ compilers
enable_language(C CXX)

set(bench_flags "-march=native -mavx") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${bench_flags}")

# rule for generating the heisenberg example program
add_executable(heisenberg heisenberg.cpp)
target_link_libraries(heisenberg ${ALPS_LIBRARIES})
 

if (ALPS_HAVE_PYTHON)

	# rule for generating python export
	set_property(GLOBAL APPEND PROPERTY PY_MODULES_LIST pyndsim)
	if(BUILD_SHARED_LIBS)
        add_library(pyndsim MODULE EXCLUDE_FROM_ALL export.cpp)
		set_target_properties(pyndsim PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
		if(WIN32 AND NOT UNIX)
			set_target_properties(pyndsim PROPERTIES SUFFIX ".pyd")
	    endif(WIN32 AND NOT UNIX)
	else(BUILD_SHARED_LIBS)
		set_property(GLOBAL APPEND PROPERTY PY_STATIC_MODULES_LIST pyndsim)
        add_library(pyndsim STATIC EXCLUDE_FROM_ALL export.cpp)
	endif (BUILD_SHARED_LIBS)
	set_target_properties(pyndsim PROPERTIES PREFIX "")
	target_link_libraries(pyndsim ${ALPS_LIBRARIES})

endif (ALPS_HAVE_PYTHON)
