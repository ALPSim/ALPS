# this is used only for building ALPS library

include $(top_builddir)config/include.mk

# overwrite variables in include.mk

install_sh = $(top_srcdir)/config/install-sh
mkinstalldirs = $(top_srcdir)/config/mkinstalldirs
ARXX = $(SHELL) $(top_builddir)config/arxx
MAKEDEPEND_PL = perl $(top_srcdir)/config/makedepend.pl

cppflags = -I$(top_builddir)src -I$(top_srcdir)/src $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
cppflags_mpi = $(cppflags) $(MPI_CPPFLAGS) -DALPS_MPI
cppflags_pvm = $(cppflags) $(PVM_CPPFLAGS) -DALPS_PVM
cxxflags = $(CXXFLAGS)
ldflags = -L$(builddir) -L$(top_builddir)src $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(XML_LDFLAGS) $(HDF5_LDFLAGS)
ldflags_mpi = $(ldflags) $(MPI_LDFLAGS)
ldflags_pvm = $(ldflags) $(PVM_LDFLAGS)
