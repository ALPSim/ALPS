DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS) -DALPS_MPI
AM_LDFLAGS = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS) $(MPI_LDFLAGS)
LDADD = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

if BUILD_UNIFIED_LIB
lib_LTLIBRARIES =
else
lib_LTLIBRARIES =
if HAVE_MPI
lib_LTLIBRARIES += libcomm-mpi.la
endif
endif

SRC = $(top_srcdir)/src/alps/osiris/buffer.C $(top_srcdir)/src/alps/osiris/dump.C $(top_srcdir)/src/alps/osiris/os.C $(top_srcdir)/src/alps/osiris/process.C $(top_srcdir)/src/alps/osiris/comm.C $(top_srcdir)/src/alps/osiris/mpdump.C

libcomm_mpi_la_SOURCES = $(SRC)
nodist_libcomm_mpi_la_SOURCES =
