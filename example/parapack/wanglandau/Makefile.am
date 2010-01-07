DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS) 

check_PROGRAMS = wanglandau
dist_noinst_DATA = params_learn params_measure params_reweight

wanglandau_SOURCES = main.C wanglandau.C wanglandau.h
wanglandau_ns_SOURCES = main_ns.C wanglandau.C wanglandau.h

TESTS = wanglandau_ns
check_PROGRAMS += $(TESTS)
dist_noinst_DATA += wanglandau_ns.ip wanglandau_ns.op
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests example examples
test tests :
example examples : check
params_learn.in.xml: params_learn
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_learn params_learn
params_measure.in.xml: params_measure
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_measure params_measure
params_reweight.in.xml: params_reweight
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_reweight params_reweight
clean-local :
	rm -f params_learn.* params_measure.* params_reweight.* weight.dump observable.dump
