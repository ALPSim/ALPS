BOOST_SRCDIR=@BOOST_SRCDIR@
DEFAULT_INCLUDES = -I$(top_builddir)/src -I$(top_srcdir)/src
AM_CPPFLAGS = $(BASE_CPPFLAGS) $(BOOST_CPPFLAGS) $(XML_CPPFLAGS) $(HDF5_CPPFLAGS) $(MPI_CPPFLAGS)

LDFLAGS_SGL = @LIB_VERSION_INFO@ $(BASE_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS) $(XML_LDFLAGS)
LIBS_SGL = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(BASE_LIBS)

LDFLAGS_MPI = $(LDFLAGS_SGL) $(MPI_LDFLAGS)
LIBS_MPI = $(BOOST_LIBS) $(HDF5_LIBS) $(XML_LIBS) $(MPI_LIBS) $(BASE_LIBS)

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
  $(top_srcdir)/src/alps/osiris/xdrcore.C \
  $(top_srcdir)/src/alps/osiris/dumparchive.C
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
src_PARAPACK_COMMON = \
  $(top_srcdir)/src/alps/parapack/clone.C \
  $(top_srcdir)/src/alps/parapack/clone_info.C \
  $(top_srcdir)/src/alps/parapack/worker_factory.C \
  $(top_srcdir)/src/alps/parapack/filelock.C \
  $(top_srcdir)/src/alps/parapack/job.C \
  $(top_srcdir)/src/alps/parapack/logger.C \
  $(top_srcdir)/src/alps/parapack/mc_worker.C \
  $(top_srcdir)/src/alps/parapack/measurement.C \
  $(top_srcdir)/src/alps/parapack/option.C \
  $(top_srcdir)/src/alps/parapack/queue.C \
  $(top_srcdir)/src/alps/parapack/rng_helper.C \
  $(top_srcdir)/src/alps/parapack/types.C \
  $(top_srcdir)/src/alps/parapack/util.C
src_PARAPACK_SGL = \
  $(top_srcdir)/src/alps/parapack/scheduler_comm.C \
  $(top_srcdir)/src/alps/parapack/scheduler_sgl.C
src_PARAPACK_MPI = \
  $(top_srcdir)/src/alps/parapack/clone_info_mpi.C \
  $(top_srcdir)/src/alps/parapack/clone_mpi.C \
  $(top_srcdir)/src/alps/parapack/parallel_factory.C \
  $(top_srcdir)/src/alps/parapack/process_mpi.C \
  $(top_srcdir)/src/alps/parapack/scheduler_comm.C \
  $(top_srcdir)/src/alps/parapack/scheduler_mpi.C

src_boost_DATE_TIME = \
  $(BOOST_SRCDIR)/date_time/src/gregorian/date_generators.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/greg_month.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/greg_weekday.cpp \
  $(BOOST_SRCDIR)/date_time/src/gregorian/gregorian_types.cpp \
  $(BOOST_SRCDIR)/date_time/src/posix_time/posix_time_types.cpp

src_boost_FILESYSTEM_1_32 = \
  $(BOOST_SRCDIR)/filesystem/src/convenience.cpp \
  $(BOOST_SRCDIR)/filesystem/src/operations_posix_windows.cpp \
  $(BOOST_SRCDIR)/filesystem/src/path_posix_windows.cpp
  # $(BOOST_SRCDIR)/filesystem/src/exception.cpp
src_boost_FILESYSTEM_1_33 = $(src_boost_FILESYSTEM_1_32)
src_boost_FILESYSTEM_1_34 = \
  $(BOOST_SRCDIR)/filesystem/src/operations.cpp \
  $(BOOST_SRCDIR)/filesystem/src/path.cpp \
  $(BOOST_SRCDIR)/filesystem/src/portability.cpp \
  $(top_builddir)/src/boost/fs_utf8_codecvt_facet.cpp
  # $(BOOST_SRCDIR)/filesystem/src/exception.cpp
src_boost_FILESYSTEM_1_35 = \
  $(BOOST_SRCDIR)/filesystem/src/operations.cpp \
  $(BOOST_SRCDIR)/filesystem/src/path.cpp \
  $(BOOST_SRCDIR)/filesystem/src/portability.cpp \
  $(top_builddir)/src/boost/fs_utf8_codecvt_facet.cpp
src_boost_FILESYSTEM_1_36 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_1_37 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_1_38 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_1_39 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_1_40 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_1_41 = $(src_boost_FILESYSTEM_1_35)
src_boost_FILESYSTEM_SVN = $(src_boost_FILESYSTEM_1_35)

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
src_boost_MPI_1_38 = $(src_boost_MPI_1_35)
src_boost_MPI_1_39 = $(src_boost_MPI_1_35)
src_boost_MPI_1_40 = $(src_boost_MPI_1_35)
src_boost_MPI_1_41 = $(src_boost_MPI_1_35)
src_boost_MPI_SVN = $(src_boost_MPI_1_35)

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
  $(BOOST_SRCDIR)/regex/src/w32_regex_traits.cpp
if BUILD_BOOST_WCHAR
src_boost_REGEX_1_33 += \
  $(BOOST_SRCDIR)/regex/src/wc_regex_traits.cpp \
  $(BOOST_SRCDIR)/regex/src/wide_posix_api.cpp \
  $(BOOST_SRCDIR)/regex/src/winstances.cpp
endif
src_boost_REGEX_1_34 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_35 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_36 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_37 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_38 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_39 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_40 = $(src_boost_REGEX_1_33)
src_boost_REGEX_1_41 = $(src_boost_REGEX_1_33)
src_boost_REGEX_SVN = $(src_boost_REGEX_1_33)

src_boost_SERIALIZATION_1_32 = \
  $(BOOST_SRCDIR)/serialization/src/basic_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_serializer_map.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_iprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_oprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_oarchive.cpp \
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
  $(BOOST_SRCDIR)/serialization/src/basic_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_archive_impl.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_serializer_map.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_iprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_oprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_xml_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/codecvt_null.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_no_rtti.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_typeid.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/stl_port.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/void_cast.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_grammar.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_oarchive.cpp
if BUILD_BOOST_WCHAR
src_boost_SERIALIZATION_1_33 += \
  $(BOOST_SRCDIR)/serialization/src/basic_text_wiprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_woprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_woarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_woarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_wgrammar.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_woarchive.cpp \
  $(top_builddir)/src/boost/sr_utf8_codecvt_facet.cpp
endif
src_boost_SERIALIZATION_1_34 = $(src_boost_SERIALIZATION_1_33)
src_boost_SERIALIZATION_1_35 = \
  $(BOOST_SRCDIR)/serialization/src/basic_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_iserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_pointer_oserializer.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_serializer_map.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_iprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_oprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_xml_archive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/codecvt_null.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_no_rtti.cpp \
  $(BOOST_SRCDIR)/serialization/src/extended_type_info_typeid.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/polymorphic_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/stl_port.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_oarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/void_cast.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_grammar.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_iarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_oarchive.cpp
if BUILD_BOOST_WCHAR
src_boost_SERIALIZATION_1_35 += \
  $(BOOST_SRCDIR)/serialization/src/basic_text_wiprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/basic_text_woprimitive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/binary_woarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/text_woarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_wgrammar.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_wiarchive.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_woarchive.cpp \
  $(top_builddir)/src/boost/sr_utf8_codecvt_facet.cpp
endif
src_boost_SERIALIZATION_1_36 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_37 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_38 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_39 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_40 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_41 = $(src_boost_SERIALIZATION_1_35)
src_boost_SERIALIZATION_1_41 = \
  $(src_boost_SERIALIZATION_1_35) \
  $(BOOST_SRCDIR)/serialization/src/archive_exception.cpp \
  $(BOOST_SRCDIR)/serialization/src/shared_ptr_helper.cpp \
  $(BOOST_SRCDIR)/serialization/src/xml_archive_exception.cpp
src_boost_SERIALIZATION_SVN = $(src_boost_SERIALIZATION_1_41)

src_boost_SIGNALS = \
  $(BOOST_SRCDIR)/signals/src/connection.cpp \
  $(BOOST_SRCDIR)/signals/src/named_slot_map.cpp \
  $(BOOST_SRCDIR)/signals/src/signal_base.cpp \
  $(BOOST_SRCDIR)/signals/src/slot.cpp \
  $(BOOST_SRCDIR)/signals/src/trackable.cpp

src_BOOST_SYSTEM_1_35 = \
  $(BOOST_SRCDIR)/system/src/error_code.cpp
src_BOOST_SYSTEM_1_36 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_1_37 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_1_38 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_1_39 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_1_40 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_1_41 = $(src_BOOST_SYSTEM_1_35)
src_BOOST_SYSTEM_SVN = $(src_BOOST_SYSTEM_1_35)

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

if BUILD_UNIFIED_LIB
lib_LTLIBRARIES = 
if ALPS_FULL
lib_LTLIBRARIES += libalps-sgl.la
if HAVE_BOOST_MPI
lib_LTLIBRARIES += libalps-mpi.la
endif
endif
else
lib_LTLIBRARIES = libalps.la
if ALPS_FULL
lib_LTLIBRARIES += libparapack-sgl.la
if HAVE_BOOST_MPI
lib_LTLIBRARIES += libparapack-mpi.la
endif
endif
endif

libalps_la_SOURCES = $(src_ALPS) $(src_EXPRESSION) $(src_OSIRIS) $(src_PARAMETER) $(src_PARSER) $(src_RANDOM) $(src_ALEA)
nodist_libalps_la_SOURCES =
libalps_la_LDFLAGS = $(LDFLAGS_SGL)
# libalps_la_LIBADD = ../comm-sgl/libcomm-sgl.la $(LIBS_SGL)

libparapack_sgl_la_SOURCES =
nodist_libparapack_sgl_la_SOURCES =
libparapack_sgl_la_LDFLAGS = $(LDFLAGS_SGL)
# libparapack_sgl_la_LIBADD = libalps.la ../comm-sgl/libcomm-sgl.la $(LIBS_SGL)

libparapack_mpi_la_SOURCES =
nodist_libparapack_mpi_la_SOURCES =
libparapack_mpi_la_LDFLAGS = $(LDFLAGS_MPI)
# libparapack_mpi_la_LIBADD = libalps.la ../comm-mpi/libcomm-mpi.la $(LIBS_MPI)

if ALPS_FULL
libalps_la_SOURCES += $(src_LATTICE) $(src_MODEL) $(src_SCHEDULER) $(src_PARAPACK_COMMON)
libparapack_sgl_la_SOURCES += $(src_PARAPACK_SGL)
if HAVE_BOOST_MPI
libparapack_mpi_la_SOURCES += $(src_PARAPACK_MPI)
endif
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

if BOOST_1_35
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_35) $(src_boost_REGEX_1_35) $(src_boost_SERIALIZATION_1_35) $(src_BOOST_SYSTEM_1_35)
endif

if BOOST_1_36
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_36) $(src_boost_REGEX_1_36) $(src_boost_SERIALIZATION_1_36) $(src_BOOST_SYSTEM_1_36)
endif

if BOOST_1_37
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_37) $(src_boost_REGEX_1_37) $(src_boost_SERIALIZATION_1_37) $(src_BOOST_SYSTEM_1_37)
endif

if BOOST_1_38
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_38) $(src_boost_REGEX_1_38) $(src_boost_SERIALIZATION_1_38) $(src_BOOST_SYSTEM_1_38)
endif

if BOOST_1_39
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_39) $(src_boost_REGEX_1_39) $(src_boost_SERIALIZATION_1_39) $(src_BOOST_SYSTEM_1_39)
endif

if BOOST_1_40
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_40) $(src_boost_REGEX_1_40) $(src_boost_SERIALIZATION_1_40) $(src_BOOST_SYSTEM_1_40)
endif

if BOOST_1_41
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_1_41) $(src_boost_REGEX_1_41) $(src_boost_SERIALIZATION_1_41) $(src_BOOST_SYSTEM_1_41)
endif

if BOOST_SVN
nodist_libalps_la_SOURCES += $(src_boost_FILESYSTEM_SVN) $(src_boost_REGEX_SVN) $(src_boost_SERIALIZATION_SVN) $(src_BOOST_SYSTEM_SVN)
endif

if BUILD_BOOST_MPI
if BOOST_1_35
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_35)
endif
if BOOST_1_36
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_36)
endif
if BOOST_1_37
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_37)
endif
if BOOST_1_38
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_38)
endif
if BOOST_1_39
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_39)
endif
if BOOST_1_40
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_40)
endif
if BOOST_1_41
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_1_41)
endif
if BOOST_SVN
nodist_libparapack_mpi_la_SOURCES += $(src_boost_MPI_SVN)
endif
endif

if BUILD_BOOST_SIGNALS
nodist_libalps_la_SOURCES += $(src_boost_SIGNALS)
endif
if BUILD_BOOST_THREAD
nodist_libalps_la_SOURCES += $(src_boost_THREAD)
endif

endif

# unified libraries

SRC_OSIRIS = $(top_srcdir)/src/alps/osiris/buffer.C $(top_srcdir)/src/alps/osiris/dump.C $(top_srcdir)/src/alps/osiris/os.C $(top_srcdir)/src/alps/osiris/process.C $(top_srcdir)/src/alps/osiris/comm.C $(top_srcdir)/src/alps/osiris/mpdump.C

libalps_sgl_la_SOURCES = $(libalps_la_SOURCES) $(libparapack_sgl_la_SOURCES) $(SRC_OSIRIS)
nodist_libalps_sgl_la_SOURCES = $(nodist_libalps_la_SOURCES) $(nodist_libparapack_sgl_la_SOURCES)
libalps_sgl_la_LDFLAGS = $(LDFLAGS_SGL)
# libalps_sgl_la_LIBADD = $(LIBS_SGL)

libalps_mpi_la_SOURCES = $(libalps_la_SOURCES) $(libparapack_mpi_la_SOURCES) $(SRC_OSIRIS)
nodist_libalps_mpi_la_SOURCES = $(nodist_libalps_la_SOURCES) $(nodist_libparapack_mpi_la_SOURCES)
libalps_mpi_la_CPPFLAGS = $(AM_CPPFLAGS) -DALPS_MPI
libalps_mpi_la_LDFLAGS = $(LDFLAGS_MPI)
# libalps_mpi_la_LIBADD = $(LIBS_MPI)

clean-local :
	rm -f *.ii *.ti
