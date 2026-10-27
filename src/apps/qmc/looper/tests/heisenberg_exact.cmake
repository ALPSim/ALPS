# The loop algorithm must reproduce the energy of small antiferromagnetic
# Heisenberg clusters that full diagonalization gives exactly, within its
# own error bars. A fixed SEED does not make the run repeat step for step,
# since the scheduler stops on elapsed time rather than on a step count, so
# this is a statistical test: five error bars.
include("${helpers}")
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

file(WRITE "${test_dir}/heisenberg" [=[
MODEL="spin"
CONSERVED_QUANTUMNUMBERS="Sz"
ALGORITHM="loop; path integral"
J=1
T=0.3
THERMALIZATION=2000
SWEEPS=20000
SEED=5
{LATTICE="dimer";}
{LATTICE="chain lattice";L=4;}
]=])

# Energy per site at T = 0.3 from fulldiag with the same MODEL, J and lattice.
set(exact_energies -0.32666221539805373 -0.47188571917054334)
set(labels "dimer" "chain L=4")

alps_run(parameter2xml "${parameter2xml}" heisenberg)
alps_run(loop "${loop}" --write-xml --Tmin 1 heisenberg.in.xml)
foreach(task RANGE 1 2)
  math(EXPR index "${task} - 1")
  list(GET labels ${index} label)
  list(GET exact_energies ${index} exact)
  alps_read_average("${test_dir}/heisenberg.task${task}.out.xml" "Energy Density")
  alps_check_within_error("Heisenberg ${label}, energy density" "${mean}" "${error}" "${exact}" 5)
endforeach()
