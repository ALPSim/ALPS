DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(SQLITE_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS) $(MPI_LDFLAGS)
if BUILD_UNIFIED_LIB
  LDADD = $(top_builddir)/lib/alps/libalps-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS)
else
  LDADD = $(top_builddir)/lib/alps/libalps.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(HDF5_LIBS) $(MPI_LIBS) $(BASE_LIBS) 
endif
BOOST_DIR = @BOOST_DIR@

check_PROGRAMS = fleas_correlated fleas_independent fleas_uncorrelated fleas_direct fleas_simpleminded

fleas_correlated_SOURCES = fleas_correlated.C fleas.h
fleas_independent_SOURCES = fleas_independent.C fleas.h
fleas_uncorrelated_SOURCES = fleas_uncorrelated.C fleas.h
fleas_direct_SOURCES = fleas_direct.C fleas.h
fleas_simpleminded_SOURCES = fleas_simpleminded.C fleas.h

.PHONY : test tests example examples
test tests :
example examples : check
