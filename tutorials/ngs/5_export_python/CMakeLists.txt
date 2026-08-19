cmake_minimum_required(VERSION 3.22)
project(alps_nanobind_export_example LANGUAGES CXX)

find_package(ALPS REQUIRED CONFIG)
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _nanobind_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
list(PREPEND CMAKE_PREFIX_PATH "${_nanobind_cmake_dir}")
find_package(nanobind 2.10 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(_alps_libraries ${ALPS_LIBRARIES})
list(TRANSFORM _alps_libraries REPLACE "^hdf5-shared$" "hdf5")
separate_arguments(_alps_link_options NATIVE_COMMAND "${ALPS_EXTRA_LINKER_FLAGS}")
separate_arguments(_alps_compile_options NATIVE_COMMAND "${ALPS_CMAKE_CXX_FLAGS}")

link_directories(${ALPS_LIBRARY_DIRS})
set(_alps_runtime_paths ${ALPS_LIBRARY_DIRS})
if(ALPS_HDF5_INCLUDE_DIR)
  get_filename_component(_hdf5_prefix "${ALPS_HDF5_INCLUDE_DIR}" DIRECTORY)
  link_directories("${_hdf5_prefix}/lib")
  list(APPEND _alps_runtime_paths "${_hdf5_prefix}/lib")
endif()

nanobind_add_module(ising_c NB_STATIC export2py.cpp ising.cpp)
target_include_directories(ising_c PRIVATE
  ${ALPS_INCLUDE_DIRS}
  ${ALPS_EXTRA_INCLUDE_DIRS})
target_compile_options(ising_c PRIVATE ${_alps_compile_options})
target_compile_definitions(ising_c PRIVATE ${ALPS_EXTRA_DEFINITIONS})
target_link_libraries(ising_c PRIVATE ${_alps_libraries})
target_link_options(ising_c PRIVATE ${_alps_link_options})
set_target_properties(ising_c PROPERTIES
  INSTALL_RPATH "${_alps_runtime_paths}"
  BUILD_RPATH "${_alps_runtime_paths}")
