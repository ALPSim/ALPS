BOOST_SRCDIR=@BOOST_SRCDIR@
DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS)
AM_LDFLAGS = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LDADD = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

src_ALPS = \
  $(top_srcdir)/src/alps/copyright.C
src_ALEA = \
  $(top_srcdir)/src/alps/alea/nan.C \
  $(top_srcdir)/src/alps/alea/observable.C \
  $(top_srcdir)/src/alps/alea/observableset.C
src_EXPRESSION = \
  $(top_srcdir)/src/alps/expression/evaluator.C 
src_OSIRIS = \
  $(top_srcdir)/src/alps/osiris/xdrdump.C \
  $(top_srcdir)/src/alps/osiris/xdrcore.C
src_LATTICE = \
  $(top_srcdir)/src/alps/lattice/disorder.C \
  $(top_srcdir)/src/alps/lattice/latticedescriptor.C \
  $(top_srcdir)/src/alps/lattice/latticegraphdescriptor.C \
  $(top_srcdir)/src/alps/lattice/latticelibrary.C \
  $(top_srcdir)/src/alps/lattice/unitcell.C
src_MODEL = \
  $(top_srcdir)/src/alps/model/modellibrary.C \
  $(top_srcdir)/src/alps/model/siteterm.C \
  $(top_srcdir)/src/alps/model/bondterm.C \
  $(top_srcdir)/src/alps/model/globaloperator.C
src_RANDOM = \
  $(top_srcdir)/src/alps/random/rngfactory.C
src_PARAMETER = \
  $(top_srcdir)/src/alps/parameter/parameter.C \
  $(top_srcdir)/src/alps/parameter/parameterlist.C \
  $(top_srcdir)/src/alps/parameter/parameters.C
src_PARSER = \
  $(top_srcdir)/src/alps/parser/parser.C \
  $(top_srcdir)/src/alps/parser/xmlattributes.C \
  $(top_srcdir)/src/alps/parser/xmlhandler.C \
  $(top_srcdir)/src/alps/parser/xmlparser.C \
  $(top_srcdir)/src/alps/parser/xmlstream.C \
  $(top_srcdir)/src/alps/parser/xslt_path.C
src_SCHEDULER = \
  $(top_srcdir)/src/alps/scheduler/abstract_task.C \
  $(top_srcdir)/src/alps/scheduler/factory.C \
  $(top_srcdir)/src/alps/scheduler/info.C \
  $(top_srcdir)/src/alps/scheduler/master_scheduler.C \
  $(top_srcdir)/src/alps/scheduler/measurement_operators.C \
  $(top_srcdir)/src/alps/scheduler/montecarlo.C \
  $(top_srcdir)/src/alps/scheduler/mpp_scheduler.C \
  $(top_srcdir)/src/alps/scheduler/options.C \
  $(top_srcdir)/src/alps/scheduler/remote_task.C \
  $(top_srcdir)/src/alps/scheduler/remote_worker.C \
  $(top_srcdir)/src/alps/scheduler/scheduler.C \
  $(top_srcdir)/src/alps/scheduler/signal.C \
  $(top_srcdir)/src/alps/scheduler/single_scheduler.C \
  $(top_srcdir)/src/alps/scheduler/serial_scheduler.C \
  $(top_srcdir)/src/alps/scheduler/slave_task.C \
  $(top_srcdir)/src/alps/scheduler/task.C \
  $(top_srcdir)/src/alps/scheduler/worker.C \
  $(top_srcdir)/src/alps/scheduler/workertask.C

src_boost_DATE_TIME = \
  $(BOOST_SRCDIR)/date_time/src/gregorian/date_generators.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/greg_month.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/greg_weekday.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/gregorian_types.cpp \
  $(BOOST_SRCDIR)/date_time/src/posix_time/posix_time_types.cpp

src_boost_FILESYSTEM_1_32 = \
  $(BOOST_SRCDIR)/filesystem/src/convenience.cpp \
  $(BOOST_SRCDIR)/filesystem/src/exception.cpp \
  $(BOOST_SRCDIR)/filesystem/src/operations_posix_windows.cpp \
  $(BOOST_SRCDIR)/filesystem/src/path_posix_windows.cpp
src_boost_FILESYSTEM_1_33 = $(src_boost_FILESYSTEM_1_32)
src_boost_FILESYSTEM_1_34 = \
  $(BOOST_SRCDIR)/filesystem/src/exception.cpp \
  $(BOOST_SRCDIR)/filesystem/src/operations.cpp \
  $(BOOST_SRCDIR)/filesystem/src/path.cpp \
  $(BOOST_SRCDIR)/filesystem/src/portability.cpp \
  $(top_builddir)/src/boost/fs_utf8_codecvt_facet.cpp
src_boost_FILESYSTEM_CVS = $(src_boost_FILESYSTEM_1_34)

src_boost_PROGRAM_OPTIONS = \
  $(BOOST_SRCDIR)/program_options/src/cmdline.cpp \
  $(BOOST_SRCDIR)/program_options/src/config_file.cpp \
  $(BOOST_SRCDIR)/program_options/src/convert.cpp \
  $(BOOST_SRCDIR)/program_options/src/options_description.cpp \
  $(BOOST_SRCDIR)/program_options/src/parsers.cpp \
  $(BOOST_SRCDIR)/program_options/src/positional_options.cpp \
  $(BOOST_SRCDIR)/program_options/src/value_semantic.cpp \
  $(BOOST_SRCDIR)/program_options/src/variables_map.cpp \
  $(BOOST_SRCDIR)/program_options/src/winmain.cpp \
  $(top_builddir)/src/boost/po_utf8_codecvt_facet.cpp

src_boost_REGEX_1_32 = \
  $(BOOST_SRCDIR)/regex/src/c_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/c_regex_traits_common.cpp \
  $(BOOST_SRCDIR)/regex/src/cpp_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/cregex.cpp \
  $(BOOST_SRCDIR)/regex/src/fileiter.cpp \
  $(BOOST_SRCDIR)/regex/src/instances.cpp \
  $(BOOST_SRCDIR)/regex/src/posix_api.cpp \
  $(BOOST_SRCDIR)/regex/src/regex.cpp \
  $(BOOST_SRCDIR)/regex/src/regex_debug.cpp \
  $(BOOST_SRCDIR)/regex/src/regex_synch.cpp \
  $(BOOST_SRCDIR)/regex/src/w32_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/wide_posix_api.cpp \
  $(BOOST_SRCDIR)/regex/src/winstances.cpp
src_boost_REGEX_1_33 = \
  $(BOOST_SRCDIR)/regex/src/c_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/cpp_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/cregex.cpp \
  $(BOOST_SRCDIR)/regex/src/fileiter.cpp \
  $(BOOST_SRCDIR)/regex/src/icu.cpp \
  $(BOOST_SRCDIR)/regex/src/instances.cpp \
  $(BOOST_SRCDIR)/regex/src/posix_api.cpp \
  $(BOOST_SRCDIR)/regex/src/regex.cpp \
  $(BOOST_SRCDIR)/regex/src/regex_debug.cpp \
  $(BOOST_SRCDIR)/regex/src/regex_raw_buffer.cpp \
  $(BOOST_SRCDIR)/regex/src/regex_traits_defaults.cpp \
  $(BOOST_SRCDIR)/regex/src/static_mutex.cpp \
  $(BOOST_SRCDIR)/regex/src/usinstances.cpp \
  $(BOOST_SRCDIR)/regex/src/w32_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/wc_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/wide_posix_api.cpp \
  $(BOOST_SRCDIR)/regex/src/winstances.cpp
src_boost_REGEX_1_34 = $(src_boost_REGEX_1_33)
src_boost_REGEX_CVS = $(src_boost_REGEX_1_33)

src_boost_SERIALIZATION_1_32 = \
  $(BOOST_SRCDIR)/serialization/src/basic_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_serializer_map.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_iprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_oprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_woarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/codecvt_null.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_no_rtti.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_typeid.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/void_cast.cpp \
  $(top_builddir)/src/boost/sr_utf8_codecvt_facet.cpp
#  $(BOOST_SRCDIR)/serialization/src/basic_text_wiprimitive.cpp
#  $(BOOST_SRCDIR)/serialization/src/basic_text_woprimitive.cpp
#  $(BOOST_SRCDIR)/serialization/src/basic_xml_archive.cpp
#  $(BOOST_SRCDIR)/serialization/src/text_wiarchive.cpp
#  $(BOOST_SRCDIR)/serialization/src/text_woarchive.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_grammar.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_iarchive.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_oarchive.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_wgrammar.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_wiarchive.cpp
#  $(BOOST_SRCDIR)/serialization/src/xml_woarchive.cpp

src_boost_SERIALIZATION_1_33 = \
  $(src_boost_SERIALIZATION_1_32) \
  $(BOOST_SRCDIR)/serialization/src/basic_archive_impl.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/stl_port.cpp
src_boost_SERIALIZATION_1_34 = $(src_boost_SERIALIZATION_1_33)
src_boost_SERIALIZATION_CVS = $(src_boost_SERIALIZATION_1_33)

src_boost_SIGNALS = \
  $(BOOST_SRCDIR)/signals/src/connection.cpp \
  $(BOOST_SRCDIR)/signals/src/named_slot_map.cpp \
  $(BOOST_SRCDIR)/signals/src/signal_base.cpp \
  $(BOOST_SRCDIR)/signals/src/slot.cpp \
  $(BOOST_SRCDIR)/signals/src/trackable.cpp

src_boost_THREAD = \
  $(BOOST_SRCDIR)/thread/src/barrier.cpp \
  $(BOOST_SRCDIR)/thread/src/condition.cpp \
  $(BOOST_SRCDIR)/thread/src/exceptions.cpp \
  $(BOOST_SRCDIR)/thread/src/mutex.cpp \
  $(BOOST_SRCDIR)/thread/src/once.cpp \
  $(BOOST_SRCDIR)/thread/src/recursive_mutex.cpp \
  $(BOOST_SRCDIR)/thread/src/thread.cpp \
  $(BOOST_SRCDIR)/thread/src/tss.cpp \
  $(BOOST_SRCDIR)/thread/src/xtime.cpp

lib_LTLIBRARIES = libalps.la 

libalps_la_SOURCES = $(src_ALPS) $(src_EXPRESSION) $(src_OSIRIS) $(src_PARAMETER) $(src_PARSER) $(src_RANDOM)
nodist_libalps_la_SOURCES =

if ALPS_FULL
libalps_la_SOURCES += $(src_ALEA) $(src_LATTICE) $(src_MODEL) $(src_SCHEDULER)
endif

if BUILD_BOOST

libalps_la_SOURCES += $(top_srcdir)/src/boost/throw_exception.C
nodist_libalps_la_SOURCES += $(src_boost_DATE_TIME) $(src_boost_PROGRAM_OPTIONS)

if BOOST_1_32
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_32) $(src_boost_REGEX_1_32) $(src_boost_SERIALIZATION_1_32)
endif

if BOOST_1_33
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_33) $(src_boost_REGEX_1_33) $(src_boost_SERIALIZATION_1_33)
endif

if BOOST_1_34
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_34) $(src_boost_REGEX_1_34) $(src_boost_SERIALIZATION_1_34)
endif

if BOOST_CVS
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_CVS) $(src_boost_REGEX_CVS) $(src_boost_SERIALIZATION_CVS)
endif

if BUILD_BOOST_SIGNALS
nodist_libalps_la_SOURCES += $(src_boost_SIGNALS)
endif
if BUILD_BOOST_THREAD
nodist_libalps_la_SOURCES += $(src_boost_THREAD)
endif

endif
