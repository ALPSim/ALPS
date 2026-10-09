# Native builds consume checked-in, validated data and need no Python parser.
set(_alps_citation_generated "${PROJECT_SOURCE_DIR}/.github/scripts/citations/generated")
include("${_alps_citation_generated}/snapshots.cmake")
if(NOT _alps_citation_generated_format EQUAL 1)
  message(FATAL_ERROR "Unsupported generated ALPS citation format")
endif()
list(LENGTH _alps_citation_files _alps_citation_count)
math(EXPR _alps_citation_last "${_alps_citation_count} - 1")
foreach(_index RANGE ${_alps_citation_last})
  list(GET _alps_citation_files ${_index} _file)
  list(GET _alps_citation_hashes ${_index} _expected)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/${_file}")
  file(SHA256 "${PROJECT_SOURCE_DIR}/${_file}" _actual)
  if(NOT _actual STREQUAL _expected)
    message(FATAL_ERROR "Generated citation data is stale (${_file}).\n"
      "After editing citation sources, run: python .github/scripts/generate_citations.py --regenerate\n"
      "Python with PyYAML/jsonschema is a maintainer tool, not a native build requirement.")
  endif()
endforeach()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_alps_citation_generated}/snapshots.cmake")

configure_file("${_alps_citation_generated}/citations_data.inc"
               "${PROJECT_BINARY_DIR}/src/alps/utility/citations_data.inc" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/CITATION.md" "${PROJECT_BINARY_DIR}/CITATION.md" COPYONLY)

install(FILES
  "${PROJECT_SOURCE_DIR}/CITATION.cff"
  "${PROJECT_SOURCE_DIR}/CITATIONS.yaml"
  "${PROJECT_BINARY_DIR}/CITATION.md"
  DESTINATION "share/alps" COMPONENT libraries)
