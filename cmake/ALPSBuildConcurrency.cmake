# SPDX-License-Identifier: MIT
# Windows targets share an output directory for their runtime DLLs.
if(WIN32 AND CMAKE_GENERATOR MATCHES "^Ninja")
  # POST_BUILD DLL staging is part of the link edge. USES_TERMINAL on that
  # command does not serialize links; give the complete edge its own pool.
  # Preserve a parent project's pools, including CMAKE_JOB_POOLS defaults.
  get_property(_alps_job_pools_set GLOBAL PROPERTY JOB_POOLS SET)
  if(NOT _alps_job_pools_set)
    set_property(GLOBAL PROPERTY JOB_POOLS "${CMAKE_JOB_POOLS}")
  endif()
  set_property(GLOBAL APPEND PROPERTY JOB_POOLS alps_runtime_link=1)
  set(CMAKE_JOB_POOL_LINK alps_runtime_link)
endif()
