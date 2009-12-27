SOURCES = abstract_fitter.C base_incr_fitter.C dummy_fitter.C est_grad_fitter.C fitting_scheduler.C spinmc_factory.C abstractspinsim.h clusterupdate.h connect.h factory.h faststack.h fitter.h fitting_scheduler.h helper.h ising.h lapack.h localupdate.h matrices.h on.h potts.h spinsim.h tinyvec.h xy.h

bin_PROGRAMS =
if BUILD_MC_SPINS
bin_PROGRAMS += spinmc spinmc_evaluate fit
endif

spinmc_SOURCES = spinmc.C $(SOURCES)
spinmc_evaluate_SOURCES = spinmc_evaluate.C $(SOURCES)
fit_SOURCES = fit.C $(SOURCES)
