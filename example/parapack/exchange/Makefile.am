DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
if BUILD_UNIFIED_LIB
  LIBS_SGL = $(top_builddir)/lib/alps/libalps-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
else
  LIBS_SGL = $(top_builddir)/lib/alps/libparapack-sgl.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
endif

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
if BUILD_UNIFIED_LIB
  LIBS_MPI = $(top_builddir)/lib/alps/libalps-mpi.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
else
  LIBS_MPI = $(top_builddir)/lib/alps/libparapack-mpi.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
endif

check_PROGRAMS = exchange
if HAVE_BOOST_MPI
check_PROGRAMS += exchange_mpi
endif
dist_noinst_DATA = params-ising params-loop

exchange_SOURCES = main.C ising.C loop.C
exchange_LDFLAGS = $(LDFLAGS_SGL)
exchange_LDADD = $(LIBS_SGL)

exchange_mpi_SOURCES = main.C ising.C loop.C ising_p.C loop_p.C
exchange_mpi_LDFLAGS = $(LDFLAGS_MPI)
exchange_mpi_LDADD = $(LIBS_MPI)

exchange_ns_SOURCES = main_ns.C ising.C loop.C
exchange_ns_LDFLAGS = $(LDFLAGS_SGL)
exchange_ns_LDADD = $(LIBS_SGL)

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
