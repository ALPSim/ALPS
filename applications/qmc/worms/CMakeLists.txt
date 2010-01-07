add_library(worm_impl STATIC WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C)
target_link_libraries(worm_impl alps boost)

add_executable(worm main.C)
add_executable(worm_evaluate evaluate.C)

target_link_libraries(worm worm_impl alps boost)
target_link_libraries(worm_evaluate worm_impl alps boost)

install(TARGETS worm RUNTIME DESTINATION bin)
install(TARGETS worm_evaluate RUNTIME DESTINATION bin)
