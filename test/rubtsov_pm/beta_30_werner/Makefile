all: test
test:
	../../../dmft < rubtsov.param
	../../difftest G_omega G_omega_converged 5.e-3
	../../difftest G_tau G_tau_converged 5.e-3
	greenclean
