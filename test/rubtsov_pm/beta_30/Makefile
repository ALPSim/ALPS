all:
	../../../dmft < rubtsov.param
	../../difftest G_omega G_omega_converged 1.e-3
	../../difftest G_tau G_tau_converged 1.e-3
	greenclean
