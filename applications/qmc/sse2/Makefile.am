SOURCES =SSE.Classes.hpp SSE.Directed.cpp SSE.Initialization.cpp SSE.Measurements.cpp SSE.Update.cpp SSE.Histogram.hpp SSE.hpp SSE.cpp

dist_noinst_DATA = conventions.txt 

bin_PROGRAMS =
if BUILD_QMC_SSE
bin_PROGRAMS += sse2 
endif

sse2_SOURCES = $(SOURCES)

