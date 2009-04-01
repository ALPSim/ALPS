DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(top_builddir)/lib/alps/libparapack-sgl.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(top_builddir)/lib/alps/libparapack-mpi.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

check_PROGRAMS = ising
if HAVE_BOOST_MPI
check_PROGRAMS += ising_mpi
endif
dist_noinst_DATA = params

ising_SOURCES = main.C ising.h ising.C
ising_LDFLAGS = $(LDFLAGS_SGL)
ising_LDADD = $(LIBS_SGL)

ising_mpi_SOURCES = main.C ising.h ising.C
ising_mpi_LDFLAGS = $(LDFLAGS_MPI)
ising_mpi_LDADD = $(LIBS_MPI)

ising_ns_SOURCES = main_ns.C ising.h ising.C
ising_ns_LDFLAGS = $(LDFLAGS_SGL)
ising_ns_LDADD = $(LIBS_SGL)

TESTS = ising_ns
check_PROGRAMS += ising_ns
dist_noinst_DATA += ising_ns.ip ising_ns.op
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

.PHONY : test tests example examples
test tests :
example examples : check
params.in.xml: $(srcdir)/params
	$(top_builddir)/tool/parameter2xml $(srcdir)/params params
clean-local :
	rm -f params.*
