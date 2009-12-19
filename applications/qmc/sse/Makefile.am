SOURCES = SSE.Classes.hpp SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp SSE.Update.cpp SSE.cpp SSE.hpp

dist_noinst_DATA = conventions.txt 

bin_PROGRAMS =
if BUILD_QMC_SSE
bin_PROGRAMS += dirloop_sse dirloop_sse_evaluate
if HAVE_MPI
bin_PROGRAMS += dirloop_sse_mpi
endif
endif

dirloop_sse_SOURCES = $(SOURCES)
dirloop_sse_evaluate_SOURCES = evaluate.C

dirloop_sse_mpi_SOURCES = $(SOURCES)
dirloop_sse_mpi_LDFLAGS = @LDFLAGS_MPI@
dirloop_sse_mpi_LDADD = @LIBS_MPI@
