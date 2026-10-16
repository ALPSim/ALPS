# Numeric locks of the classic DMRG component against ALPS::ed on small
# systems, where enough kept states make the sweeps exact.
foreach(name IN ITEMS energies measurements validation operator_forms scratch)
  add_executable(classic_dmrg_${name} ${name}.cpp)
  target_link_libraries(classic_dmrg_${name} PRIVATE ALPS::classic_dmrg ALPS::ed)
  alps_add_gtest(classic_dmrg_${name} COMPONENT classic_dmrg LABELS "classic_dmrg;numeric" TIMEOUT 600)
endforeach()
