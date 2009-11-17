#
# libboost
#

set(BOOST_SOURCES "")
set(BOOST_SHARED_CPPFLAGS ${BOOST_SHARED_CPPFLAGS} BOOST_ALL_DYN_LINK=1)
set(BOOST_STATIC_CPPFLAGS ${BOOST_STATIC_CPPFLAGS} BOOST_ALL_STATIC_LINK=1)

# Boost.Date_Time
if(BUILD_BOOST_DATE_TIME)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/date_time/src")
  set(SOURCES gregorian/date_generators.cpp gregorian/greg_month.cpp gregorian/greg_weekday.cpp
      gregorian/gregorian_types.cpp posix_time/posix_time_types.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_DATE_TIME)

# Boost.Filesystem
if(BUILD_BOOST_FILESYSTEM)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/filesystem/src")
  set(SOURCES
    operations.cpp path.cpp portability.cpp utf8_codecvt_facet.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_FILESYSTEM)

# Boost.ProgramOptions
if(BUILD_BOOST_PROGRAM_OPTIONS)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/program_options/src")
  set(SOURCES cmdline.cpp config_file.cpp options_description.cpp parsers.cpp
      variables_map.cpp value_semantic.cpp positional_options.cpp utf8_codecvt_facet.cpp
      convert.cpp winmain.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_PROGRAM_OPTIONS)

# Boost.Regex
if(BUILD_BOOST_REGEX)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/regex/src")
  set(SOURCES c_regex_traits.cpp cpp_regex_traits.cpp cregex.cpp fileiter.cpp icu.cpp
      instances.cpp posix_api.cpp regex.cpp regex_debug.cpp regex_raw_buffer.cpp
      regex_traits_defaults.cpp static_mutex.cpp usinstances.cpp w32_regex_traits.cpp
      wc_regex_traits.cpp wide_posix_api.cpp winstances.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_REGEX)

# Boost.Serialization
if(BUILD_BOOST_SERIALIZATION)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/serialization/src")
  set(SOURCES basic_archive.cpp basic_iarchive.cpp basic_iserializer.cpp basic_oarchive.cpp
      basic_oserializer.cpp basic_pointer_iserializer.cpp basic_pointer_oserializer.cpp
      basic_serializer_map.cpp basic_text_iprimitive.cpp basic_text_oprimitive.cpp
      basic_text_wiprimitive.cpp basic_text_woprimitive.cpp basic_xml_archive.cpp
      binary_iarchive.cpp binary_oarchive.cpp binary_wiarchive.cpp
      binary_woarchive.cpp codecvt_null.cpp extended_type_info.cpp
      extended_type_info_no_rtti.cpp extended_type_info_typeid.cpp polymorphic_iarchive.cpp
      polymorphic_oarchive.cpp stl_port.cpp text_iarchive.cpp text_oarchive.cpp text_wiarchive.cpp
      text_woarchive.cpp utf8_codecvt_facet.cpp void_cast.cpp xml_grammar.cpp xml_iarchive.cpp
      xml_oarchive.cpp xml_wgrammar.cpp xml_wiarchive.cpp xml_woarchive.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_SERIALIZATION)

# Boost.System
if(BUILD_BOOST_SYSTEM)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/system/src")
  set(SOURCES error_code.cpp)
  foreach(S ${SOURCES})
    set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_SYSTEM)

add_library(boost ${BOOST_SOURCES})
if(BUILD_SHARED_LIBS)
  set_target_properties(boost PROPERTIES COMPILE_DEFINITIONS "${BOOST_SHARED_CPPFLAGS}")
else(BUILD_SHARED_LIBS)
  set_target_properties(boost PROPERTIES COMPILE_DEFINITIONS "${BOOST_STATIC_CPPFLAGS}")
endif(BUILD_SHARED_LIBS)

#######################################################################
# install
#######################################################################
install(TARGETS boost
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
