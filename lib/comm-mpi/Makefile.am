DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS) -DALPS_MPI
AM_LDFLAGS = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS) $(MPI_LDFLAGS)
LDADD = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

lib_LTLIBRARIES =
if HAVE_MPI
lib_LTLIBRARIES += libcomm-mpi.la
endif

SRC = $(top_srcdir)/src/alps/osiris/buffer.C $(top_srcdir)/src/alps/osiris/dump.C $(top_srcdir)/src/alps/osiris/os.C $(top_srcdir)/src/alps/osiris/process.C $(top_srcdir)/src/alps/osiris/comm.C $(top_srcdir)/src/alps/osiris/mpdump.C
if HAVE_BOOST_MPI
SRC += $(top_srcdir)/src/alps/parapack/clone_info_mpi.C \
  $(top_srcdir)/src/alps/parapack/clone_mpi.C \
  $(top_srcdir)/src/alps/parapack/parallel_factory.C \
  $(top_srcdir)/src/alps/parapack/process_mpi.C \
  $(top_srcdir)/src/alps/parapack/scheduler_mpi.C
endif

src_boost_MPI_1_35 = \
  $(BOOST_SRCDIR)/mpi/src/broadcast.cpp \
  $(BOOST_SRCDIR)/mpi/src/communicator.cpp \
  $(BOOST_SRCDIR)/mpi/src/computation_tree.cpp \
  $(BOOST_SRCDIR)/mpi/src/content_oarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/environment.cpp \
  $(BOOST_SRCDIR)/mpi/src/exception.cpp \
  $(BOOST_SRCDIR)/mpi/src/graph_communicator.cpp \
  $(BOOST_SRCDIR)/mpi/src/group.cpp \
  $(BOOST_SRCDIR)/mpi/src/intercommunicator.cpp \
  $(BOOST_SRCDIR)/mpi/src/mpi_datatype_cache.cpp \
  $(BOOST_SRCDIR)/mpi/src/mpi_datatype_oarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/packed_iarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/packed_oarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/packed_skeleton_iarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/packed_skeleton_oarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/point_to_point.cpp \
  $(BOOST_SRCDIR)/mpi/src/request.cpp \
  $(BOOST_SRCDIR)/mpi/src/text_skeleton_oarchive.cpp \
  $(BOOST_SRCDIR)/mpi/src/timer.cpp
# $(BOOST_SRCDIR)/mpi/src/text_skeleton_iarchive.cpp
src_boost_MPI_1_36 = $(src_boost_MPI_1_35)
src_boost_MPI_1_37 = $(src_boost_MPI_1_35)
src_boost_MPI_SVN = $(src_boost_MPI_1_35)

libcomm_mpi_la_SOURCES = $(SRC)
nodist_libcomm_mpi_la_SOURCES =
if BUILD_BOOST_MPI
if BOOST_1_35
nodist_libcomm_mpi_la_SOURCES += $(src_boost_MPI_1_35)
endif
if BOOST_1_36
nodist_libcomm_mpi_la_SOURCES += $(src_boost_MPI_1_36)
endif
if BOOST_1_37
nodist_libcomm_mpi_la_SOURCES += $(src_boost_MPI_1_37)
endif
if BOOST_SVN
nodist_libcomm_mpi_la_SOURCES += $(src_boost_MPI_SVN)
endif
endif
