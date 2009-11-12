include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tool)
include_directories(${Boost_ROOT_DIR})

add_executable(ising main.C ising.C)
add_dependencies(ising alps)
target_link_libraries(ising alps)
install(TARGETS ising
        RUNTIME DESTINATION bin)

