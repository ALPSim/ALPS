DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(SQLITE_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS)
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

check_PROGRAMS = wanglandau
if HAVE_BOOST_MPI
check_PROGRAMS += wanglandau_mpi
endif
dist_noinst_DATA = params_learn params_measure params_reweight

wanglandau_SOURCES = main.C wanglandau.C wanglandau.h
wanglandau_LDFLAGS = $(LDFLAGS_SGL)
wanglandau_LDADD = $(LIBS_SGL)

wanglandau_mpi_SOURCES = main.C wanglandau.C wanglandau.h
wanglandau_mpi_LDFLAGS = $(LDFLAGS_MPI)
wanglandau_mpi_LDADD = $(LIBS_MPI)

wanglandau_ns_SOURCES = main_ns.C wanglandau.C wanglandau.h
wanglandau_ns_LDFLAGS = $(LDFLAGS_SGL)
wanglandau_ns_LDADD = $(LIBS_SGL)

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
