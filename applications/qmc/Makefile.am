SUBDIRS =
if EXIST_QMC_CHECKSIGN
SUBDIRS += checksign
endif
if EXIST_QMC_LOOPER
SUBDIRS += looper
endif
if EXIST_QMC_SSE
SUBDIRS += sse
endif
if EXIST_QMC_SSE_LOOP
SUBDIRS += sse-loop
endif
if EXIST_QMC_WORMS
SUBDIRS += worms
endif
if EXIST_QMC_QWL
SUBDIRS += qwl
endif
if EXIST_QMC_SSE2
SUBDIRS += sse2
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

noinst_HEADERS = qmc.h
