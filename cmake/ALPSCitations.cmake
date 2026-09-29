# Validate the citation authorities and embed an offline snapshot in libalps.
# Python is a build tool here; Python bindings/development headers are not needed.
if(NOT ALPS_CITATION_PYTHON)
  if(PYTHON_INTERPRETER)
    set(_alps_citation_python "${PYTHON_INTERPRETER}")
  else()
    find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
    set(_alps_citation_python "${Python3_EXECUTABLE}")
  endif()
  set(ALPS_CITATION_PYTHON "${_alps_citation_python}" CACHE FILEPATH
      "Python interpreter with PyYAML and jsonschema for citation generation")
endif()

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/CITATION.cff"
  "${PROJECT_SOURCE_DIR}/CITATIONS.yaml"
  "${PROJECT_SOURCE_DIR}/script/generate_citations.py"
  "${PROJECT_SOURCE_DIR}/script/citations/cff-1.2.0.schema.json"
  "${PROJECT_SOURCE_DIR}/script/citations/policy.schema.json")

execute_process(
  COMMAND "${ALPS_CITATION_PYTHON}" "${PROJECT_SOURCE_DIR}/script/generate_citations.py"
          --root "${PROJECT_SOURCE_DIR}"
          --cpp "${PROJECT_BINARY_DIR}/src/alps/utility/citations_data.inc"
          --markdown "${PROJECT_BINARY_DIR}/CITATION.md"
          --snapshot-cpp "${PROJECT_BINARY_DIR}/src/alps/utility/citation_snapshots.inc"
          --software-version "${ALPS_VERSION}"
  RESULT_VARIABLE _alps_citation_result
  ERROR_VARIABLE _alps_citation_error)
if(NOT _alps_citation_result EQUAL 0)
  message(FATAL_ERROR "Cannot generate ALPS citations using ${ALPS_CITATION_PYTHON}:\n"
    "${_alps_citation_error}")
endif()

if(ALPS_PYTHON_WHEEL)
  set(_alps_citation_destination "pyalps/share/alps")
else()
  set(_alps_citation_destination "share/alps")
endif()
install(FILES
  "${PROJECT_SOURCE_DIR}/CITATION.cff"
  "${PROJECT_SOURCE_DIR}/CITATIONS.yaml"
  "${PROJECT_BINARY_DIR}/CITATION.md"
  DESTINATION "${_alps_citation_destination}" COMPONENT libraries)
