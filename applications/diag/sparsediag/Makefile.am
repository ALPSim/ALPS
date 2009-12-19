SOURCES = sparsediag.C sparsediag.h factory.C factory.h ../matrix.h

dist_noinst_DATA =

bin_PROGRAMS =
if BUILD_DIAG_SPARSEDIAG
bin_PROGRAMS += sparsediag
if HAVE_MPI
bin_PROGRAMS += sparsediag_mpi
endif
endif

sparsediag_SOURCES = $(SOURCES)

sparsediag_mpi_SOURCES = $(SOURCES)
sparsediag_mpi_LDFLAGS = @LDFLAGS_MPI@
sparsediag_mpi_LDADD = @LIBS_MPI@
