DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS) 

check_PROGRAMS =
if HAVE_BOOST_MPI
check_PROGRAMS += ising ising_ns
endif
dist_noinst_DATA = params

ising_SOURCES = main.C ising.h ising.C
ising_ns_SOURCES = main_ns.C ising.h ising.C

TESTS =
check_PROGRAMS +=
dist_noinst_DATA += ising_ns.ip ising_ns.op-2 ising_ns.op-4
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests example examples
test tests :
example examples : check
params.in.xml: $(srcdir)/params
	$(top_builddir)/tool/parameter2xml $(srcdir)/params params
clean-local :
	rm -f params.*
