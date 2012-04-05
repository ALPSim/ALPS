include_directories(${PROJECT_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${Boost_ROOT_DIR})


FOREACH (name   example_accumulator example_custom_accum example_histogram)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} alps)
    add_alps_test(${name})
ENDFOREACH(name)
