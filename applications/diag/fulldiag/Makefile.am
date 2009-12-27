SOURCES = factory.C factory.h fulldiag.C fulldiag.h
SOURCES_EVAL = fulldiag.h fulldiag_evaluate.C

dist_noinst_DATA =

bin_PROGRAMS =
if BUILD_DIAG_FULLDIAG
bin_PROGRAMS += fulldiag fulldiag_evaluate
endif

fulldiag_SOURCES = $(SOURCES)
fulldiag_evaluate_SOURCES =$(SOURCES_EVAL)
