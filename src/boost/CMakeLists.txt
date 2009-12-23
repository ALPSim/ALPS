#
# libboost
#

set(BOOST_SOURCES "")

# Boost.Date_Time
if(BUILD_BOOST_DATE_TIME)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/date_time/src")
  set(SOURCES gregorian/date_generators.cpp gregorian/greg_month.cpp gregorian/greg_weekday.cpp
      gregorian/gregorian_types.cpp posix_time/posix_time_types.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_DATE_TIME)

# Boost.Filesystem
if(BUILD_BOOST_FILESYSTEM)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/filesystem/src")
  set(SOURCES
    operations.cpp path.cpp portability.cpp utf8_codecvt_facet.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_FILESYSTEM)

# Boost.ProgramOptions
if(BUILD_BOOST_PROGRAM_OPTIONS)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/program_options/src")
  set(SOURCES cmdline.cpp config_file.cpp options_description.cpp parsers.cpp
      variables_map.cpp value_semantic.cpp positional_options.cpp utf8_codecvt_facet.cpp
      convert.cpp winmain.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
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
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_REGEX)

# Boost.Serialization
if(BUILD_BOOST_SERIALIZATION)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/serialization/src")
  set(SOURCES archive_exception.cpp basic_archive.cpp
      basic_iarchive.cpp basic_iserializer.cpp basic_oarchive.cpp
      basic_oserializer.cpp basic_pointer_iserializer.cpp
      basic_pointer_oserializer.cpp basic_serializer_map.cpp
      basic_text_iprimitive.cpp basic_text_oprimitive.cpp
      basic_text_wiprimitive.cpp basic_text_woprimitive.cpp
      basic_xml_archive.cpp binary_iarchive.cpp binary_oarchive.cpp
      binary_wiarchive.cpp binary_woarchive.cpp codecvt_null.cpp
      extended_type_info.cpp extended_type_info_no_rtti.cpp
      extended_type_info_typeid.cpp polymorphic_iarchive.cpp
      polymorphic_oarchive.cpp shared_ptr_helper.cpp stl_port.cpp
      text_iarchive.cpp text_oarchive.cpp text_wiarchive.cpp
      text_woarchive.cpp utf8_codecvt_facet.cpp void_cast.cpp
      xml_archive_exception.cpp xml_grammar.cpp xml_iarchive.cpp
      xml_oarchive.cpp xml_wgrammar.cpp xml_wiarchive.cpp
      xml_woarchive.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_SERIALIZATION)

# Boost.Python
if(BUILD_BOOST_PYTHON)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/python/src")
  set(SOURCES dict.cpp errors.cpp exec.cpp import.cpp list.cpp long.cpp
      module.cpp numeric.cpp object_operators.cpp object_protocol.cpp slice.cpp
      str.cpp tuple.cpp wrapper.cpp converter/arg_to_python_base.cpp
      converter/builtin_converters.cpp converter/from_python.cpp
      converter/registry.cpp converter/type_id.cpp object/class.cpp
      object/enum.cpp object/function.cpp object/function_doc_signature.cpp
      object/inheritance.cpp object/iterator.cpp object/life_support.cpp
      object/pickle_support.cpp object/stl_iterator.cpp
  )
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_PYTHON)

# Boost.System
if(BUILD_BOOST_SYSTEM)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/system/src")
  set(SOURCES error_code.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_SYSTEM)


add_library(boost ${BOOST_SOURCES})

if (ALPS_FOR_VISTRAILS)
  SET_TARGET_PROPERTIES(boost
    PROPERTIES BUILD_WITH_INSTALL_RPATH ON
    INSTALL_NAME_DIR "@executable_path/../Frameworks"
  )
endif (ALPS_FOR_VISTRAILS)



#######################################################################
# install
#######################################################################
if (ALPS_FOR_VISTRAILS)
install(TARGETS boost
          RUNTIME DESTINATION .
          ARCHIVE DESTINATION .
          LIBRARY DESTINATION .)
else (ALPS_FOR_VISTRAILS)
  install(TARGETS boost
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
endif (ALPS_FOR_VISTRAILS)
