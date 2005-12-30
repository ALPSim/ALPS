DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src

CPPFLAGS_SGL = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

CPPFLAGS_MPI = $(CPPFLAGS_SGL) $(MPI_CPPFLAGS)
LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

TESTS = 
TESTS_norun = evaluate evaluate2 main main2
if HAVE_MPI
TESTS_norun += main_mpi main2_mpi
endif
check_PROGRAMS = $(TESTS) $(TESTS_norun)
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)

evaluate_SOURCES = evaluate.C
evaluate_CPPFLAGS = $(CPPFLAGS_SGL)
evaluate_LDFLAGS = $(LDFLAGS_SGL)
evaluate_LDADD = $(LIBS_SGL)

evaluate2_SOURCES = evaluate2.C
evaluate2_CPPFLAGS = $(CPPFLAGS_SGL)
evaluate2_LDFLAGS = $(LDFLAGS_SGL)
evaluate2_LDADD = $(LIBS_SGL)

main_SOURCES = main.C ising.h ising.C
main_CPPFLAGS = $(CPPFLAGS_SGL)
main_LDFLAGS = $(LDFLAGS_SGL)
main_LDADD = $(LIBS_SGL)

main2_SOURCES = main2.C ising2.h ising2.C
main2_CPPFLAGS = $(CPPFLAGS_SGL)
main2_LDFLAGS = $(LDFLAGS_SGL)
main2_LDADD = $(LIBS_SGL)

main_mpi_SOURCES = main.C ising.h ising.C
main_mpi_CPPFLAGS = $(CPPFLAGS_MPI)
main_mpi_LDFLAGS = $(LDFLAGS_MPI)
main_mpi_LDADD = $(LIBS_MPI)

main2_mpi_SOURCES = main2.C ising2.h ising2.C
main2_mpi_CPPFLAGS = $(CPPFLAGS_MPI)
main2_mpi_LDFLAGS = $(LDFLAGS_MPI)
main2_mpi_LDADD = $(LIBS_MPI)

dist_noinst_DATA = parm parm2

.PHONY : test tests example examples
test tests :
example examples : check
