DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

lib_LTLIBRARIES = libcomm-sgl.la

SRC = $(top_srcdir)/src/alps/osiris/buffer.C $(top_srcdir)/src/alps/osiris/dump.C $(top_srcdir)/src/alps/osiris/os.C $(top_srcdir)/src/alps/osiris/process.C $(top_srcdir)/src/alps/osiris/comm.C $(top_srcdir)/src/alps/osiris/mpdump.C \
  $(top_srcdir)/src/alps/parapack/scheduler_sgl.C

libcomm_sgl_la_SOURCES = $(SRC)
