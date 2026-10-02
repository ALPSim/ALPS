cmake_minimum_required(VERSION 3.24)
project(ALPSReconciliationProbe LANGUAGES C CXX)

# Configure twice, in separate build directories and with separate SDK prefixes.
# These providers have overlapping headers and symbols: never link them together.
set(RECONCILIATION_PROVIDER "ALPS" CACHE STRING "ALPS or ALPSCore")
set_property(CACHE RECONCILIATION_PROVIDER PROPERTY STRINGS ALPS ALPSCore)
add_executable(reconciliation_probe probe.cpp)
target_compile_features(reconciliation_probe PRIVATE cxx_std_17)
if(RECONCILIATION_PROVIDER STREQUAL "ALPS")
  find_package(ALPS CONFIG REQUIRED)
  target_link_libraries(reconciliation_probe PRIVATE ALPS::params)
elseif(RECONCILIATION_PROVIDER STREQUAL "ALPSCore")
  find_package(ALPSCore CONFIG REQUIRED COMPONENTS params hdf5)
  target_compile_definitions(reconciliation_probe PRIVATE PROBE_ALPSCORE)
  target_link_libraries(reconciliation_probe PRIVATE alps::alps-params)
else()
  message(FATAL_ERROR "RECONCILIATION_PROVIDER must be ALPS or ALPSCore")
endif()
