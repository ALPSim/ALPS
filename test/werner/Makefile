all: test
test:
	../../dmft < werner.param
	../difftest G_omega G_omega_converged 5.e-4
	../difftest G_tau G_tau_converged 5.e-4
	greenclean
