DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(SQLITE_CPPFLAGS) $(XML_CPPFLAGS) 
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(SQLITE_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(SQLITE_LIBS) $(XML_LIBS) $(BASE_LIBS)
BOOST_DIR = @BOOST_DIR@

pkgdatadir = $(libdir)/xml

bin_PROGRAMS = compactrun convert2xml parameter2xml
dist_bin_SCRIPTS = archivecat xslttransform transformall 
bin_SCRIPTS = convert2text convert2html extractgp extracthtml extracttext extractxmgr plot2gp plot2html plot2text plot2xmgr changestylesheet

compactrun_SOURCES = compactrun.C
convert2xml_SOURCES = convert2xml.C
parameter2xml_SOURCES = parameter2xml.C

if BUILD_ARCHIVE
bin_PROGRAMS += archive
endif
archive_SOURCES = archive.cpp archive_index.cpp archive_index.hpp archive_node.cpp archive_node.hpp archive_plot.cpp archive_plot.hpp archive_sqlite.cpp archive_sqlite.hpp archive_xml.cpp archive_xml.hpp

if BUILD_DOC
bin_PROGRAMS += quickbook
nodist_quickbook_SOURCES = $(BOOST_DIR)/tools/quickbook/detail/quickbook.cpp 
if BOOST_1_34
nodist_quickbook_SOURCES += \
 $(BOOST_DIR)/tools/quickbook/detail/actions.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/post_process.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/utils.cpp
endif
if BOOST_1_35
nodist_quickbook_SOURCES += \
 $(BOOST_DIR)/tools/quickbook/detail/actions.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/actions_class.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/collector.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/post_process.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/template_stack.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/utils.cpp
endif
if BOOST_SVN
nodist_quickbook_SOURCES += \
 $(BOOST_DIR)/tools/quickbook/detail/actions.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/actions_class.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/collector.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/post_process.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/template_stack.cpp \
 $(BOOST_DIR)/tools/quickbook/detail/utils.cpp
endif
endif

