DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

pkgdatadir = $(libdir)/xml

bin_PROGRAMS = compactrun convert2xml parameter2xml
dist_bin_SCRIPTS = archivecat xslttransform transformall 
bin_SCRIPTS = convert2text convert2html extractgp extracthtml extracttext extractxmgr plot2gp plot2html plot2text plot2xmgr changestylesheet

compactrun_SOURCES = compactrun.C
convert2xml_SOURCES = convert2xml.C
parameter2xml_SOURCES = parameter2xml.C
