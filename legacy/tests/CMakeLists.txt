# Build-tree checks for the deprecated C++ interfaces. Installed-SDK consumers
# (UseALPS.cmake, include.mk) are exercised by tests/cmake/test_sdk_contracts.py.
add_executable(legacy_params_path_constructor params_path_constructor.cpp)
target_link_libraries(legacy_params_path_constructor PRIVATE ALPS::alps)
alps_add_test(legacy_params_path_constructor)

add_executable(legacy_save_observable_header save_observable_header.cpp)
target_link_libraries(legacy_save_observable_header PRIVATE ALPS::alps)
alps_add_test(legacy_save_observable_header)

set_tests_properties(legacy_params_path_constructor legacy_save_observable_header
  PROPERTIES LABELS "legacy")
