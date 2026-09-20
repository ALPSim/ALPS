cmake_minimum_required(VERSION 3.27)
project(pyalps_parameter_contracts LANGUAGES C CXX)
find_package(ALPS REQUIRED CONFIG)
find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _nanobind_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
list(PREPEND CMAKE_PREFIX_PATH "${_nanobind_dir}")
find_package(nanobind 2.10 CONFIG REQUIRED)
set(CMAKE_CXX_STANDARD 17)
nanobind_add_module(parameter_probe NB_STATIC probe.cpp)
alps_target_link_pyalps(parameter_probe PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
