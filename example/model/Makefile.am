DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(top_builddir)/lib/alps/libalps.la $(top_builddir)/lib/comm-sgl/libcomm-sgl.la $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

check_PROGRAMS = print_numeric print_symbolic

print_numeric_SOURCES = print_numeric.C matrix.h 
print_symbolic_SOURCES = print_symbolic.C matrix.h 

noinst_HEADERS = matrix.h
dist_noinst_DATA = parm_numeric parm_symbolic

.PHONY : test tests example examples
test tests :
example examples : check
