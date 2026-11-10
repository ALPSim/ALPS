# Adapters cross the typed params and legacy parameter/parser libraries.
add_executable(param_legacy_adapters legacy_adapters.cpp)
target_link_libraries(param_legacy_adapters PRIVATE ALPS::alps)
alps_add_gtest(param_legacy_adapters COMPONENT params LABELS integration compatibility)

add_executable(param_parameter_input parameter_input.cpp)
target_link_libraries(param_parameter_input PRIVATE ALPS::alps)
alps_add_gtest(param_parameter_input COMPONENT params LABELS integration compatibility)
