SOURCES = externalsolver.C hirschfyesim.C fouriertransform.C hilberttransformer.C alps_solver.C 2dsimpson.C selfconsistency.C auxiliaryfunctions.C 
INCLUDES = -I@top_srcdir@/dmft/
SOURCES_HF = solver_main.C
dist_noinst_DATA =

bin_PROGRAMS =
bin_PROGRAMS += dmft hirschfye
if HAVE_MPI
bin_PROGRAMS += dmft_mpi
endif
dmft_SOURCES = $(SOURCES) main.C
dmft_mpi_SOURCES = $(SOURCES) main.C
dmft_mpi_LDFLAGS = @LDFLAGS_MPI@
dmft_mpi_LDADD = @LIBS_MPI@
hirschfye_SOURCES= $(SOURCES_HF) $(SOURCES)

if BUILD_DMFT_SOLVERS_HYBRIDIZATION
INCLUDES +=-DBUILD_DMFT_SOLVERS_HYBRIDIZATION
dmft_SOURCES += ../DMFT-solvers/Werner/impurity.C
dmft_SOURCES += ../DMFT-solvers/Werner/fourpoint.C
dmft_mpi_SOURCES += ../DMFT-solvers/Werner/impurity.C
dmft_mpi_SOURCES += ../DMFT-solvers/Werner/fourpoint.C
endif
