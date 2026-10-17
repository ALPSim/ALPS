# Numeric locks of the classic DMRG component against ALPS::ed on small
# systems, where enough kept states make the sweeps exact.
foreach(name IN ITEMS energies measurements validation operator_forms scratch mps)
  add_executable(classic_dmrg_${name} ${name}.cpp)
  target_link_libraries(classic_dmrg_${name} PRIVATE ALPS::classic_dmrg ALPS::ed)
  alps_add_gtest(classic_dmrg_${name} COMPONENT classic_dmrg LABELS "classic_dmrg;numeric" TIMEOUT 600)
endforeach()

# Time evolution against dense Trotter products and exact evolution on the
# ED sector basis (Eigen's dense solver).
add_executable(classic_dmrg_evolution evolution.cpp)
target_link_libraries(classic_dmrg_evolution PRIVATE ALPS::classic_dmrg ALPS::ed Eigen3::Eigen)
alps_add_gtest(classic_dmrg_evolution COMPONENT classic_dmrg LABELS "classic_dmrg;numeric" TIMEOUT 600)

# Unevenly distributed ground states from a warm-up configuration against
# sorted fields and free fermions (Eigen's dense solver).
add_executable(classic_dmrg_warmup warmup.cpp)
target_link_libraries(classic_dmrg_warmup PRIVATE ALPS::classic_dmrg ALPS::ed Eigen3::Eigen)
alps_add_gtest(classic_dmrg_warmup COMPONENT classic_dmrg LABELS "classic_dmrg;numeric" TIMEOUT 600)

# Classic molecular DMRG against an occupation-bit exact diagonalization
# (Eigen's dense solver) and against ALPS::ed.
add_executable(classic_dmrg_molecular molecular.cpp)
target_link_libraries(classic_dmrg_molecular PRIVATE ALPS::classic_dmrg ALPS::ed Eigen3::Eigen)
alps_add_gtest(classic_dmrg_molecular COMPONENT classic_dmrg LABELS "classic_dmrg;numeric;chemistry" TIMEOUT 600)

# The private engine on four-site fermion products. It compiles the
# header-only engine itself, so it needs the private sources, Eigen and BLAS.
add_executable(classic_dmrg_multisite multisite.cpp)
target_include_directories(classic_dmrg_multisite PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../src")
target_link_libraries(classic_dmrg_multisite PRIVATE Eigen3::Eigen ${ALPS_BLAS_TARGET})
alps_add_gtest(classic_dmrg_multisite COMPONENT classic_dmrg LABELS "classic_dmrg;numeric" TIMEOUT 600)
