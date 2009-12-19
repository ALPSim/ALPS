SOURCES = loop.C sse.C sse.h vertex.C vertex.h

dist_noinst_DATA = 

noinst_PROGRAMS =
if BUILD_QMC_SSE_LOOP
noinst_PROGRAMS += loop
if HAVE_MPI
noinst_PROGRAMS += loop_mpi
endif
endif

loop_SOURCES = $(SOURCES)

loop_mpi_SOURCES = $(SOURCES)
loop_mpi_LDFLAGS = @LDFLAGS_MPI@
loop_mpi_LDADD = @LIBS_MPI@
