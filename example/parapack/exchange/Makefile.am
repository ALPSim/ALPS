DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS) 

check_PROGRAMS = exchange
dist_noinst_DATA = params-ising params-loop

exchange_SOURCES = main.C ising.C loop.C
exchange_ns_SOURCES = main_ns.C ising.C loop.C

TESTS = exchange_ns
check_PROGRAMS += $(TESTS)
dist_noinst_DATA += exchange_ns.ip exchange_ns.op
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests example examples
test tests :
example examples : check
params-ising.in.xml:
	$(top_builddir)/tool/parameter2xml $(srcdir)/params-ising params-ising
params-loop.in.xml:
	$(top_builddir)/tool/parameter2xml $(srcdir)/params-loop params-loop
clean-local :
	rm -f params-ising.* params-loop.*
