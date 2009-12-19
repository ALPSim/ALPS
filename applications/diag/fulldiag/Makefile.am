SOURCES = factory.C factory.h fulldiag.C fulldiag.h
SOURCES_EVAL = fulldiag.h fulldiag_evaluate.C

dist_noinst_DATA =

bin_PROGRAMS =
if BUILD_DIAG_FULLDIAG
bin_PROGRAMS += fulldiag fulldiag_evaluate
if HAVE_MPI
bin_PROGRAMS += fulldiag_mpi
endif
endif

fulldiag_SOURCES = $(SOURCES)
fulldiag_evaluate_SOURCES =$(SOURCES_EVAL)

fulldiag_mpi_SOURCES = $(SOURCES)
fulldiag_mpi_LDFLAGS = @LDFLAGS_MPI@
fulldiag_mpi_LDADD = @LIBS_MPI@
