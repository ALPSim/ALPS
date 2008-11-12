DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(top_builddir)/lib/alps/libparapack-sgl.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(top_builddir)/lib/alps/libparapack-mpi.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

check_PROGRAMS = main_ps
if HAVE_BOOST_MPI
check_PROGRAMS += main_pp
endif
dist_noinst_DATA = params_learn params_measure params_reweight

main_ps_SOURCES = main.C wanglandau.C wanglandau.h
main_ps_LDFLAGS = $(LDFLAGS_SGL)
main_ps_LDADD = $(LIBS_SGL)

main_pp_SOURCES = main.C wanglandau.C wanglandau.h
main_pp_LDFLAGS = $(LDFLAGS_MPI)
main_pp_LDADD = $(LIBS_MPI)

main_ns_SOURCES = main_ns.C wanglandau.C wanglandau.h
main_ns_LDFLAGS = $(LDFLAGS_SGL)
main_ns_LDADD = $(LIBS_SGL)

TESTS = main_ns
check_PROGRAMS += $(TESTS)
dist_noinst_DATA += main_ns.ip main_ns.op
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)
.PHONY : test tests
test tests : check

.PHONY : params_learn params_measure params_reweight
params : params_learn.in.xml params_measure.in.xml params_reweight.in.xml
params_learn.in.xml: params_learn
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_learn params_learn
params_measure.in.xml: params_measure
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_measure params_measure
params_reweight.in.xml: params_reweight
	$(top_builddir)/tool/parameter2xml $(srcdir)/params_reweight params_reweight
clean-local :
	rm -f params_learn.* params_measure.* params_reweight.* weight.dump observable.dump
