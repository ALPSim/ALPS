cmake_minimum_required(VERSION 3.27)
project(alps_nanobind_export_example LANGUAGES C CXX)

find_package(ALPS REQUIRED CONFIG)
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _nanobind_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
list(PREPEND CMAKE_PREFIX_PATH "${_nanobind_cmake_dir}")
find_package(nanobind 2.10 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

nanobind_add_module(ising_c NB_STATIC export2py.cpp ising.cpp)
alps_target_link_pyalps(ising_c PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
