SOURCES =SSE.Classes.hpp SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp SSE.Update.cpp SSE.Histogram.hpp SSE.hpp SSE.cpp

dist_noinst_DATA = conventions.txt 

bin_PROGRAMS =
if BUILD_QMC_SSE
bin_PROGRAMS += sse2 
if HAVE_MPI
bin_PROGRAMS += sse2_mpi 
endif
endif

sse2_SOURCES = $(SOURCES)
sse2_mpi_SOURCES = $(SOURCES)
sse2_mpi_LDFLAGS = @LDFLAGS_MPI@
sse2_mpi_LDADD = @LIBS_MPI@

