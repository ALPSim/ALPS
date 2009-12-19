dist_noinst_DATA = 
SUBDIRS =
if EXIST_TEST_FULLDIAG
SUBDIRS += fulldiag
endif
if EXIST_TEST_SPARSEDIAG
SUBDIRS += sparsediag
endif
if EXIST_TEST_CHECKSIGN
SUBDIRS += checksign
endif
if EXIST_TEST_LOOPER
SUBDIRS += looper
endif
if EXIST_TEST_SSE
SUBDIRS += sse
endif
if EXIST_TEST_SSE_LOOP
SUBDIRS += sse-loop
endif
if EXIST_TEST_WORMS
SUBDIRS += worms
endif
if EXIST_TEST_QWL
SUBDIRS += qwl
endif
if EXIST_TEST_MC
SUBDIRS += mc
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

.PHONY : test tests example examples
test tests example examples : check
