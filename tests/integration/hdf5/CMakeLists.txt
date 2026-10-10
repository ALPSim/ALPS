# Cross-component archive compatibility with legacy parameters and observables.
foreach(name hdf5_parms hdf5_ising hdf5_observableset)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_gtest(${name} COMPONENT hdf5 LABELS integration compatibility TIMEOUT 120)
endforeach()
