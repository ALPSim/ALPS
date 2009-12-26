add_executable(qwl qwl.C)
add_executable(qwl_evaluate qwl_evaluate.C)

target_link_libraries(qwl alps boost)
target_link_libraries(qwl_evaluate alps boost)
install(TARGETS qwl RUNTIME DESTINATION bin)
install(TARGETS qwl_evaluate RUNTIME DESTINATION bin)
