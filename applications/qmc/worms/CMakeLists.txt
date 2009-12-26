add_executable(worm main.C
               WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C
)
add_executable(worm_evaluate evaluate.C
               WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C
)

target_link_libraries(worm alps boost)
target_link_libraries(worm_evaluate alps boost)
install(TARGETS worm RUNTIME DESTINATION bin)
install(TARGETS worm_evaluate RUNTIME DESTINATION bin)
