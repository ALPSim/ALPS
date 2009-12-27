DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(SQLITE_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)

if BUILD_UNIFIED_LIB
  LDADD = $(top_builddir)/lib/alps/libalps-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS)
else
  LDADD = $(top_builddir)/lib/alps/libalps.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS) 
endif

check_PROGRAMS = evaluate evaluate2 main main2 main3

evaluate_SOURCES = evaluate.C
evaluate2_SOURCES = evaluate2.C
main_SOURCES = main.C ising.h ising.C
main2_SOURCES = main2.C ising2.h ising2.C
main3_SOURCES = main3.C ising.h ising.C

dist_noinst_DATA = parm parm2 parm3

.PHONY : test tests example examples
test tests :
example examples : check
