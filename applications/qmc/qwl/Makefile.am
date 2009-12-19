SOURCES = qwl.C qwl_sse.h qwl_histogram.h
SOURCES_EVAL = qwl_evaluate.C

dist_noinst_DATA = 

bin_PROGRAMS =
if BUILD_QMC_QWL
bin_PROGRAMS += qwl qwl_evaluate
if HAVE_MPI
bin_PROGRAMS += qwl_mpi
endif
endif

qwl_SOURCES = $(SOURCES)

qwl_evaluate_SOURCES =$(SOURCES_EVAL)

qwl_mpi_SOURCES = $(SOURCES)
qwl_mpi_LDFLAGS = @LDFLAGS_MPI@
qwl_mpi_LDADD = @LIBS_MPI@
