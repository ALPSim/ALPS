SOURCES = WModel.C Wcheck.C Winit.C WRun.C Wmeas.C Wdostep.C main.C WKink.h WRun.h cyclic_iterator.h random.h time_struct.h

dist_noinst_DATA = 

bin_PROGRAMS =
if BUILD_QMC_WORMS
bin_PROGRAMS += worm worm_evaluate
endif

worm_SOURCES = $(SOURCES)
worm_evaluate_SOURCES = evaluate.C

