DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
if BUILD_UNIFIED_LIB
  LIBS_SGL = $(top_builddir)/lib/alps/libalps-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
else
  LIBS_SGL = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
endif

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
if BUILD_UNIFIED_LIB
  LIBS_MPI = $(top_builddir)/lib/alps/libalps-mpi.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
else
  LIBS_MPI = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
endif

check_PROGRAMS = evaluate evaluate2 main main2 main3
if HAVE_MPI
check_PROGRAMS += main_mpi main2_mpi main3_mpi
endif

evaluate_SOURCES = evaluate.C
evaluate_LDFLAGS = $(LDFLAGS_SGL)
evaluate_LDADD = $(LIBS_SGL)

evaluate2_SOURCES = evaluate2.C
evaluate2_LDFLAGS = $(LDFLAGS_SGL)
evaluate2_LDADD = $(LIBS_SGL)

main_SOURCES = main.C ising.h ising.C
main_LDFLAGS = $(LDFLAGS_SGL)
main_LDADD = $(LIBS_SGL)

main2_SOURCES = main2.C ising2.h ising2.C
main2_LDFLAGS = $(LDFLAGS_SGL)
main2_LDADD = $(LIBS_SGL)

main3_SOURCES = main3.C ising.h ising.C
main3_LDFLAGS = $(LDFLAGS_SGL)
main3_LDADD = $(LIBS_SGL)

main_mpi_SOURCES = main.C ising.h ising.C
main_mpi_LDFLAGS = $(LDFLAGS_MPI)
main_mpi_LDADD = $(LIBS_MPI)

main2_mpi_SOURCES = main2.C ising2.h ising2.C
main2_mpi_LDFLAGS = $(LDFLAGS_MPI)
main2_mpi_LDADD = $(LIBS_MPI)

main3_mpi_SOURCES = main3.C ising.h ising.C
main3_mpi_LDFLAGS = $(LDFLAGS_MPI)
main3_mpi_LDADD = $(LIBS_MPI)

dist_noinst_DATA = parm parm2 parm3

.PHONY : test tests example examples
test tests :
example examples : check
