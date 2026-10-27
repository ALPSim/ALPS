# The two-dimensional Ising model on small periodic square lattices has an
# exactly known energy. spinmc must reproduce it, with both update schemes,
# within its own error bars. A fixed SEED does not make the run repeat
# step for step, since the scheduler stops on elapsed time rather than on
# a step count, so this is a statistical test: five error bars.
include("${helpers}")
file(REMOVE_RECURSE "${test_dir}")
file(MAKE_DIRECTORY "${test_dir}")
set(ENV{ALPS_XML_PATH} "${xml_dir}")

# Total energy <E> of the L x L periodic Ising model, H = -J sum s_i s_j,
# at T = 2.269186 (the critical temperature of the infinite lattice), J = 1,
# from the transfer matrix.
set(sizes 2 4 8)
set(exact_energies -6.788224 -25.049971 -95.461647)

foreach(update cluster local)
  set(name "ising_${update}")
  file(WRITE "${test_dir}/${name}" "LATTICE=\"square lattice\"
MODEL=\"Ising\"
J=1
T=2.269186
UPDATE=\"${update}\"
THERMALIZATION=10000
SWEEPS=50000
SEED=10
{L=2;}
{L=4;}
{L=8;}
")
  alps_run(parameter2xml "${parameter2xml}" "${name}")
  alps_run(spinmc "${spinmc}" --Tmin 1 "${name}.in.xml")
  foreach(task RANGE 1 3)
    math(EXPR index "${task} - 1")
    list(GET sizes ${index} L)
    list(GET exact_energies ${index} exact)
    alps_run(convert2xml "${convert2xml}" "${name}.task${task}.out.xml")
    alps_read_average("${test_dir}/${name}.task${task}.out.xml" "Energy")
    alps_check_within_error("Ising L=${L}, ${update} updates, energy" "${mean}" "${error}" "${exact}" 5)
  endforeach()
endforeach()
