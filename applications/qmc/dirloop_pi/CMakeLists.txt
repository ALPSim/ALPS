add_executable(dirloop_pi mail.cpp lowa.cpp)
target_link_libraries(dirloop_pi alps boost)
install(TARGETS dirloop_pi RUNTIME DESTINATION bin COMPONENT applications)
