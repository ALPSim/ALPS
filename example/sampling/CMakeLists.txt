#
# examples of MC error analysis
# 

if (ALPS_BUILD_EXAMPLES)
  set (EXAMPLES fleas_correlated fleas_direct fleas_independent 
       fleas_simpleminded fleas_uncorrelated)
     
  foreach(program ${EXAMPLES})
    add_executable(${program} ${program}.C)
    target_link_libraries(${program} alps)
    install(TARGETS ${program} COMPONENT examples RUNTIME DESTINATION bin )
  endforeach(program)
endif (ALPS_BUILD_EXAMPLES)

