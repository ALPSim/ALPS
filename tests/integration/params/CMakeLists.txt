# Adapters cross the typed params and legacy parameter/parser libraries.
add_executable(param_legacy_adapters legacy_adapters.cpp)
target_link_libraries(param_legacy_adapters PRIVATE ALPS::alps)
alps_add_gtest(param_legacy_adapters COMPONENT params LABELS integration compatibility)

# A path must not reach the HDF5 constructor through implicit conversion.
add_executable(param_path_overload path_overload.cpp)
target_link_libraries(param_path_overload PRIVATE ALPS::alps)
alps_add_gtest(param_path_overload COMPONENT params LABELS integration compatibility)
