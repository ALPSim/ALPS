# Run alps2fcidump on INPUT and compare its FCIDUMP with REFERENCE byte for byte.
# The test couplings are dyadic, so every integral is exact in binary.
execute_process(COMMAND "${TOOL}" "${INPUT}" out.fcidump RESULT_VARIABLE status)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "alps2fcidump failed: ${status}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files out.fcidump "${REFERENCE}"
  RESULT_VARIABLE differs)
if(differs)
  file(READ out.fcidump written)
  message(FATAL_ERROR "out.fcidump differs from ${REFERENCE}:\n${written}")
endif()
