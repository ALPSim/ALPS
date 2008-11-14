DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(top_builddir)/lib/alps/libparapack-sgl.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(top_builddir)/lib/alps/libparapack-mpi.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

check_PROGRAMS = main_ps
if HAVE_BOOST_MPI
check_PROGRAMS += main_pp
endif
dist_noinst_DATA = params

main_ps_SOURCES = main.C ising.h ising.C
main_ps_LDFLAGS = $(LDFLAGS_SGL)
main_ps_LDADD = $(LIBS_SGL)

main_pp_SOURCES = main.C ising.h ising.C
main_pp_LDFLAGS = $(LDFLAGS_MPI)
main_pp_LDADD = $(LIBS_MPI)

main_ns_SOURCES = main_ns.C ising.h ising.C
main_ns_LDFLAGS = $(LDFLAGS_SGL)
main_ns_LDADD = $(LIBS_SGL)

TESTS = main_ns
check_PROGRAMS += main_ns
dist_noinst_DATA += main_ns.ip main_ns.op
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests example examples
test tests :
example examples : check
params.in.xml: $(srcdir)/params
	$(top_builddir)/tool/parameter2xml $(srcdir)/params params
clean-local :
	rm -f params.*
