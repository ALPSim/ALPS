add_executable(application_input application_input.cpp)
target_link_libraries(application_input PRIVATE ALPS::alps)
alps_add_gtest(application_input COMPONENT params)
