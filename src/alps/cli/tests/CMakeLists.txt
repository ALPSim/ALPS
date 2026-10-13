add_executable(cli_contract cli_contract.cpp)
target_link_libraries(cli_contract PRIVATE ALPS::cli)
alps_add_test(cli_contract)
set_property(TEST cli_contract PROPERTY LABELS cli)
