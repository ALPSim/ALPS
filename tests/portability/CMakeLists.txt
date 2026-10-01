add_executable(portability portability.cpp)
target_link_libraries(portability PRIVATE ALPS::alps)
add_test(NAME portability COMMAND portability)
add_executable(numerics numerics.cpp)
target_link_libraries(numerics PRIVATE ALPS::alps)
add_test(NAME numerics COMMAND numerics)

if(ALPS_ENABLE_OPENMP)
  add_executable(openmp_worker openmp_worker.cpp)
  target_link_libraries(openmp_worker PRIVATE ALPS::alps)
  add_test(NAME openmp_worker COMMAND openmp_worker)
  set_tests_properties(openmp_worker PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_DYNAMIC=FALSE;OMP_THREAD_LIMIT=2"
    TIMEOUT 30)
endif()
