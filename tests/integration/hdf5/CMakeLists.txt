# Cross-component archive compatibility with legacy parameters and observables.
foreach(name hdf5_parms hdf5_ising hdf5_observableset)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS "hdf5;integration")
endforeach()
