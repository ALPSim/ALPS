SOURCES = WModel.C Wcheck.C Winit.C WRun.C Wmeas.C Wdostep.C main.C WKink.h WRun.h cyclic_iterator.h random.h time_struct.h

dist_noinst_DATA = 

bin_PROGRAMS =
if BUILD_QMC_WORMS
bin_PROGRAMS += worm worm_evaluate
if HAVE_MPI
bin_PROGRAMS += worm_mpi
endif
endif

worm_SOURCES = $(SOURCES)

worm_evaluate_SOURCES = evaluate.C

worm_mpi_SOURCES = $(SOURCES)
worm_mpi_LDFLAGS = @LDFLAGS_MPI@
worm_mpi_LDADD = @LIBS_MPI@
