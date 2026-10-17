# Time evolution (adaptive time-dependent DMRG)

Real- and imaginary-time evolution with `pyalps.classic_dmrg.evolve`:
quenches, a melting domain wall, ramped drives and Loschmidt echoes. The
initial state of every run is the DMRG ground state of an initial
Hamiltonian; it evolves under a final Hamiltonian by second-order Trotter
gates applied while sweeping (White and Feiguin, Phys. Rev. Lett. 93,
076401 (2004)). `method="targeting"` evolves by time-step targeting
instead (Feiguin and White, Phys. Rev. B 72, 020404(R) (2005)), for final
Hamiltonians beyond nearest neighbours or larger time steps.

Read the lessons in this order:

1. [Time evolution basics (imaginary and real time)](time_evolution_intro.py)
2. [Interaction quench of a hard-core boson chain (Loschmidt echo)](interaction_quench.py)
3. [Domain-wall ("kink") dynamics of an XX spin-1/2 chain](domain_wall.py)
4. [Sudden and ramped quenches of the hard-core boson chain](quench_hardcore_boson.py)

Each script runs in a few seconds, checks its physics with assertions and
ends with `All physics assertions passed.`; `pyalps/tests/test_time_evolution_tutorials.py`
runs them. They are jupytext percent scripts, so `jupytext --to notebook`
turns them into notebooks.

These lessons are ports of the time-evolution tutorials of
skilledwolf/alpsx (`docs/tutorials/time_evolution/` at d37e51227), which
used TEBD on matrix product states; two of those were in turn ports of the
legacy ALPS TEBD tutorials (`bhquench`, `kink`). Each script's header says
what changed.
