BOOST_SRCDIR=@BOOST_SRCDIR@
DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

src_ALPS = $(top_srcdir)/src/alps/evaluator.C $(top_srcdir)/src/alps/expression.C $(top_srcdir)/src/alps/parameterlist.C $(top_srcdir)/src/alps/parameters.C $(top_srcdir)/src/alps/copyright.C
src_ALEA = $(top_srcdir)/src/alps/alea/nan.C $(top_srcdir)/src/alps/alea/observable.C $(top_srcdir)/src/alps/alea/observableset.C
src_OSIRIS = $(top_srcdir)/src/alps/osiris/dump.C $(top_srcdir)/src/alps/osiris/xdrdump.C
src_LATTICE = $(top_srcdir)/src/alps/lattice/disorder.C $(top_srcdir)/src/alps/lattice/latticedescriptor.C $(top_srcdir)/src/alps/lattice/latticegraphdescriptor.C $(top_srcdir)/src/alps/lattice/latticelibrary.C $(top_srcdir)/src/alps/lattice/unitcell.C
src_MODEL = $(top_srcdir)/src/alps/model/modellibrary.C
src_RANDOM = $(top_srcdir)/src/alps/random/rngfactory.C
src_PARSER = $(top_srcdir)/src/alps/parser/parser.C $(top_srcdir)/src/alps/parser/xmlstream.C $(top_srcdir)/src/alps/parser/xmlhandler.C $(top_srcdir)/src/alps/parser/xmlparser.C $(top_srcdir)/src/alps/parser/xslt.C
src_SCHEDULER = $(top_srcdir)/src/alps/scheduler/abstract_task.C $(top_srcdir)/src/alps/scheduler/remote_task.C $(top_srcdir)/src/alps/scheduler/remote_worker.C $(top_srcdir)/src/alps/scheduler/scheduler.C $(top_srcdir)/src/alps/scheduler/info.C $(top_srcdir)/src/alps/scheduler/single_scheduler.C $(top_srcdir)/src/alps/scheduler/master_scheduler.C $(top_srcdir)/src/alps/scheduler/slave_task.C $(top_srcdir)/src/alps/scheduler/montecarlo.C $(top_srcdir)/src/alps/scheduler/task.C $(top_srcdir)/src/alps/scheduler/mpp_scheduler.C $(top_srcdir)/src/alps/scheduler/worker.C $(top_srcdir)/src/alps/scheduler/signal.C $(top_srcdir)/src/alps/scheduler/options.C $(top_srcdir)/src/alps/scheduler/workertask.C $(top_srcdir)/src/alps/scheduler/factory.C

src_boost_DATE_TIME = $(BOOST_SRCDIR)/date_time/src/gregorian/date_generators.cpp $(BOOST_SRCDIR)/date_time/src/gregorian/greg_month.cpp $(BOOST_SRCDIR)/date_time/src/gregorian/greg_weekday.cpp $(BOOST_SRCDIR)/date_time/src/gregorian/gregorian_types.cpp $(BOOST_SRCDIR)/date_time/src/posix_time/posix_time_types.cpp
src_boost_FILESYSTEM = $(BOOST_SRCDIR)/filesystem/src/convenience.cpp $(BOOST_SRCDIR)/filesystem/src/exception.cpp $(BOOST_SRCDIR)/filesystem/src/operations_posix_windows.cpp $(BOOST_SRCDIR)/filesystem/src/path_posix_windows.cpp
src_boost_REGEX = $(BOOST_SRCDIR)/regex/src/c_regex_traits.cpp $(BOOST_SRCDIR)/regex/src/c_regex_traits_common.cpp $(BOOST_SRCDIR)/regex/src/cpp_regex_traits.cpp $(BOOST_SRCDIR)/regex/src/cregex.cpp $(BOOST_SRCDIR)/regex/src/fileiter.cpp $(BOOST_SRCDIR)/regex/src/instances.cpp $(BOOST_SRCDIR)/regex/src/posix_api.cpp $(BOOST_SRCDIR)/regex/src/regex.cpp $(BOOST_SRCDIR)/regex/src/regex_debug.cpp $(BOOST_SRCDIR)/regex/src/regex_synch.cpp $(BOOST_SRCDIR)/regex/src/w32_regex_traits.cpp $(BOOST_SRCDIR)/regex/src/wide_posix_api.cpp $(BOOST_SRCDIR)/regex/src/winstances.cpp
src_boost_SIGNALS = \
  $(BOOST_SRCDIR)/signals/src/connection.cpp \
  $(BOOST_SRCDIR)/signals/src/signal_base.cpp \
  $(BOOST_SRCDIR)/signals/src/slot.cpp \
  $(BOOST_SRCDIR)/signals/src/trackable.cpp
src_boost_SIGNALS_1_32 = \
  $(BOOST_SRCDIR)/signals/src/named_slot_map.cpp
src_boost_THREAD = \
  $(BOOST_SRCDIR)/thread/src/condition.cpp \
  $(BOOST_SRCDIR)/thread/src/exceptions.cpp \
  $(BOOST_SRCDIR)/thread/src/mutex.cpp \
  $(BOOST_SRCDIR)/thread/src/once.cpp \
  $(BOOST_SRCDIR)/thread/src/recursive_mutex.cpp \
  $(BOOST_SRCDIR)/thread/src/thread.cpp \
  $(BOOST_SRCDIR)/thread/src/threadmon.cpp \
  $(BOOST_SRCDIR)/thread/src/tss.cpp \
  $(BOOST_SRCDIR)/thread/src/xtime.cpp
src_boost_THREAD_1_32 = \
  $(BOOST_SRCDIR)/thread/src/barrier.cpp

lib_LTLIBRARIES = libalps.la

libalps_la_SOURCES = $(src_ALPS) $(src_OSIRIS) $(src_PARSER) $(src_RANDOM)
nodist_libalps_la_SOURCES =
if ALPS_FULL
libalps_la_SOURCES += $(src_ALEA) $(src_LATTICE) $(src_MODEL) $(src_SCHEDULER)
endif

if BUILD_BOOST
libalps_la_SOURCES += $(top_srcdir)/src/boost/throw_exception.C
nodist_libalps_la_SOURCES += $(src_boost_DATE_TIME) $(src_boost_FILESYSTEM) $(src_boost_REGEX)
endif
if BUILD_BOOST_SIGNALS
nodist_libalps_la_SOURCES += $(src_boost_SIGNALS)
if BOOST_1_32
nodist_libalps_la_SOURCES += $(src_boost_SIGNALS_1_32)
endif
endif
if BUILD_BOOST_THREAD
nodist_libalps_la_SOURCES += $(src_boost_THREAD)
if BOOST_1_32
nodist_libalps_la_SOURCES += $(src_boost_THREAD_1_32)
endif
endif
