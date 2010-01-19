add_executable(checksign checksign.C)
target_link_libraries(checksign alps boost)
install(TARGETS checksign RUNTIME DESTINATION bin COMPONENT applications)
