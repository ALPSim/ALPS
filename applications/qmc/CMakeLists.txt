add_subdirectory(looper)
add_subdirectory(sse)
add_subdirectory(checksign)
add_subdirectory(worms)
add_subdirectory(qwl)

IF(EXISTS ${PROJECT_SOURCE_DIR}/applications/qmc/sse-loop)
  add_subdirectory(sse-loop)
ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/applications/qmc/sse-loop)
