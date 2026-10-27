# Quantum Wang-Landau sampling gives the energy of the Heisenberg chain over
# a whole temperature range from one run. It must agree with full
# diagonalization to a few percent. Wang-Landau results carry no error bar;
# with 20 refinement steps, seeds 1 to 9 all land within 1.6% of the exact
# energies at the three temperatures checked, so the tolerance is 3%.
include("${helpers}")
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

file(WRITE "${test_dir}/heisenberg" [=[
MODEL="spin"
LATTICE="chain lattice"
L=10
J=1
T_MIN=0.5
T_MAX=1.5
DELTA_T=0.5
CUTOFF=50
NUMBER_OF_WANG_LANDAU_STEPS=20
SEED=4
]=])

# Energy per site of the L=10 chain from fulldiag at T = 0.5, 1.0, 1.5.
set(exact_energies -0.341775 -0.204655 -0.137616)

alps_run(parameter2xml "${parameter2xml}" heisenberg)
alps_run(qwl "${qwl}" --Tmin 1 heisenberg.in.xml)
alps_run(qwl_evaluate "${qwl_evaluate}" heisenberg.task1.out.xml)
alps_read_plot("${test_dir}/heisenberg.task1.plot.energy.xml")
list(LENGTH ys count)
if(NOT count EQUAL 3)
  message(FATAL_ERROR "expected 3 temperatures in the energy plot, got ${count}: ${xs}")
endif()
foreach(index RANGE 0 2)
  list(GET xs ${index} T)
  list(GET ys ${index} energy)
  list(GET exact_energies ${index} exact)
  alps_check_within_percent("Heisenberg chain L=10, T=${T}, energy density" "${energy}" "${exact}" 3)
endforeach()
