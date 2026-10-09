add_executable(cli_contract cli_contract.cpp)
target_link_libraries(cli_contract PRIVATE ALPS::cli)
alps_add_gtest(cli_contract COMPONENT cli LABELS compatibility)
