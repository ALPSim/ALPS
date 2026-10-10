# Use the historical stdin and golden-output fixtures without a test-framework
# migration. Each launch has its own working directory for lock/dump files.
file(MAKE_DIRECTORY "${OUTPUT}")
set(input "${SOURCE}/${NAME}.ip-${RANKS}")
set(expected "${SOURCE}/${NAME}.op-${RANKS}")
set(input_args "")
if(EXISTS "${input}")
  list(APPEND input_args INPUT_FILE "${input}")
endif()
execute_process(COMMAND "${MPIEXEC}" "${NUMPROC_FLAG}" "${RANKS}" ${PREFLAGS}
  "${EXECUTABLE}" ${POSTFLAGS} ${input_args}
  WORKING_DIRECTORY "${OUTPUT}" OUTPUT_FILE "${OUTPUT}/stdout.actual"
  ERROR_VARIABLE error RESULT_VARIABLE result TIMEOUT 280)
if(NOT result STREQUAL "0")
  message(FATAL_ERROR "${NAME}/${RANKS} failed (${result}): ${error}; see ${OUTPUT}/stdout.actual")
endif()
if(EXISTS "${expected}")
  execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files --ignore-eol
    "${expected}" "${OUTPUT}/stdout.actual" RESULT_VARIABLE mismatch)
  if(mismatch)
    message(FATAL_ERROR "${NAME}/${RANKS}: ${OUTPUT}/stdout.actual differs from ${expected}")
  endif()
endif()
