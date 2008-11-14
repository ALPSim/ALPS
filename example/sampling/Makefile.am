DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)

LDFLAGS_SGL = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-mpi/libcomm-mpi.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

check_PROGRAMS = fleas_correlated fleas_independent fleas_uncorrelated fleas_direct fleas_simpleminded

fleas_correlated_SOURCES = fleas_correlated.C fleas.h
fleas_correlated_LDFLAGS = $(LDFLAGS_SGL)
fleas_correlated_LDADD = $(LIBS_SGL)

fleas_independent_SOURCES = fleas_independent.C fleas.h
fleas_independent_LDFLAGS = $(LDFLAGS_SGL)
fleas_independent_LDADD = $(LIBS_SGL)

fleas_uncorrelated_SOURCES = fleas_uncorrelated.C fleas.h
fleas_uncorrelated_LDFLAGS = $(LDFLAGS_SGL)
fleas_uncorrelated_LDADD = $(LIBS_SGL)

fleas_direct_SOURCES = fleas_direct.C fleas.h
fleas_direct_LDFLAGS = $(LDFLAGS_SGL)
fleas_direct_LDADD = $(LIBS_SGL)

fleas_simpleminded_SOURCES = fleas_simpleminded.C fleas.h
fleas_simpleminded_LDFLAGS = $(LDFLAGS_SGL)
fleas_simpleminded_LDADD = $(LIBS_SGL)

.PHONY : test tests example examples
test tests :
example examples : check
