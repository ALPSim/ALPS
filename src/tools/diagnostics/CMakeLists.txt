add_executable(pconfig pconfig.C)
target_link_libraries(pconfig PRIVATE ALPS::utilities)
install(TARGETS pconfig
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

alps_register_module(NAME tools_diagnostics
  TARGETS pconfig
  DEPENDS configuration utilities)
