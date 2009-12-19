dist_noinst_DATA = 
SUBDIRS =
if EXIST_TEST_MC_SPINS
SUBDIRS += spins
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

.PHONY : test tests example examples
test tests example examples : check
