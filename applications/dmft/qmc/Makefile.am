SOURCES = externalsolver.C hirschfyesim.C fouriertransform.C hilberttransformer.C alps_solver.C 2dsimpson.C selfconsistency.C auxiliaryfunctions.C 
INCLUDES = -I.
SOURCES_HF = solver_main.C
dist_noinst_DATA =

bin_PROGRAMS =
bin_PROGRAMS += dmft hirschfye
dmft_SOURCES = $(SOURCES) main.C
hirschfye_SOURCES= $(SOURCES_HF) $(SOURCES)

if BUILD_DMFT_QMC_HYBRIDIZATION
INCLUDES +=-DBUILD_DMFT_QMC_HYBRIDIZATION
dmft_SOURCES += hybridization/impurity.C
dmft_SOURCES += hybridization/fourpoint.C
endif
