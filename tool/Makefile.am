DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS)
if BUILD_UNIFIED_LIB
  LDADD = $(top_builddir)/lib/alps/libalps-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
else
  LDADD = $(top_builddir)/lib/alps/libparapack-sgl.la $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS) $(HDF5_LIBS)
endif
BOOST_DIR = @BOOST_DIR@

pkgdatadir = $(libdir)/xml

bin_PROGRAMS = compactrun convert2xml parameter2xml printgraph pconfig pevaluate poutput txt2archive xml2archive
dist_bin_SCRIPTS = archivecat xslttransform transformall use_local_stylesheet
bin_SCRIPTS = convert2text convert2html extractgp extracthtml extracttext extractxmgr plot2gp plot2html plot2text plot2xmgr changestylesheet extractmpl plot2mpl
nobase_pkgpython_PYTHON =
dist_noinst_DATA = lattice2xml.ip lattice2xml.op
if HAVE_VTKPYTHON
  bin_PROGRAMS += lattice2xml
  bin_SCRIPTS += lattice-preview
  nobase_pkgpython_PYTHON += config.py license.py preview.py
endif

compactrun_SOURCES = compactrun.C
convert2xml_SOURCES = convert2xml.C
lattice2xml_SOURCES = lattice2xml.C
parameter2xml_SOURCES = parameter2xml.C
printgraph_SOURCES = printgraph.C
pconfig_SOURCES = pconfig.C
pevaluate_SOURCES = pevaluate.C
poutput_SOURCES = poutput.C
txt2archive_SOURCES = txt2archive.C
xml2archive_SOURCES = xml2archive.C

if BUILD_ARCHIVE
bin_PROGRAMS += archive
endif
archive_SOURCES = archive.cpp archive_index.cpp archive_index.hpp archive_node.cpp archive_node.hpp archive_plot.cpp archive_plot.hpp archive_sqlite.cpp archive_sqlite.hpp archive_xml.cpp archive_xml.hpp

if BUILD_DOC
bin_PROGRAMS += quickbook
nodist_quickbook_SOURCES = $(BOOST_DIR)/tools/quickbook/detail/quickbook.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/actions.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/actions_class.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/collector.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/input_path.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/markups.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/post_process.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/template_stack.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/utils.cpp
endif

TESTS =
if HAVE_VTKPYTHON
  TESTS += lattice2xml
endif
TESTS_ENVIRONMENT = $(top_srcdir)/config/run-test $(srcdir)
