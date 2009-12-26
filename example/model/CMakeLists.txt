#
# examples of using the model library for numeric and symbolic computations
# 

if (ALPS_BUILD_EXAMPLES)
  set (EXAMPLES print_numeric print_symbolic)
  set (INPUT parm_numeric parm_symbolic)
  
  foreach(file ${INPUT})
    configure_file(${file}.in ${CMAKE_CURRENT_BINARY_DIR}/${file})
  endforeach(file)
     
  foreach(program ${EXAMPLES})
    add_executable(${program} ${program}.C)
    target_link_libraries(${program} alps)
  endforeach(program)
endif (ALPS_BUILD_EXAMPLES)

