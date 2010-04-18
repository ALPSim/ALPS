# - Find LPSolve
# Find LPSolve headers and libraries.
#
#  LPSolve_FOUND        - True if LPSolve found
#  LPSolve_INCLUDE_DIR  - where to find lpkit.h
#  LPSolve_LIBRARIES    - list of libraries when using LPSolve
#  LPSolve_DLLS         - list of DLLs when using LPSolve

if(NOT LPSolve_INCLUDE_DIR)
  set(__TRIAL_PATHS "$ENV{LPSOLVE_ROOT}/include" "${LPSOLVE_ROOT}" "${LPSOLVE_ROOT}/include"
      /usr/include /usr/local/include /opt/include /usr/local/lpsolve /sw/include /opt/local/include/lpsolve
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/include" "$ENV{HOME}/src/lp_solve_4.0"
     )
  find_path(LPSolve_INCLUDE_DIR lpkit.h ${__TRIAL_PATHS})
  set(__TRIAL_PATHS "$ENV{LPSOLVE_ROOT}/lib" "${LPSOLVE_ROOT}" "${LPSOLVE_ROOT}/lib"
      /usr/lib /usr/local/lib /opt/lib /opt/local/lib /usr/local/lpsolve /sw/lib
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/lib" "$ENV{HOME}/src/lp_solve_4.0"
     )
  find_library(LPSolve_LIBRARIES lpsolve ${__TRIAL_PATHS})
  if(NOT LPSolve_LIBRARIES)
    find_library(LPSolve_LIBRARIES lpsolve55 ${__TRIAL_PATHS})
  endif(NOT LPSolve_LIBRARIES)
  if(NOT LPSolve_LIBRARIES)
    find_library(LPSolve_LIBRARIES lpk5 ${__TRIAL_PATHS})
  endif(NOT LPSolve_LIBRARIES)
  if(NOT LPSolve_LIBRARIES)
    find_library(LPSolve_LIBRARIES lpk ${__TRIAL_PATHS})
  endif(NOT LPSolve_LIBRARIES)
  set(__TRIAL_PATHS "$ENV{LPSOLVE_ROOT}/bin" "${LPSOLVE_ROOT}" "${LPSOLVE_ROOT}/bin"
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/bin"
     )
  find_path(LPSolve_DLL_DIR lpsolve.dll ${__TRIAL_PATHS})
  if(LPSolve_DLL_DIR)
    set(LPSolve_DLLS "${LPSolve_DLL_DIR}/lpsolve.dll" CACHE PATH "LPSolve DLL")
  endif(LPSolve_DLL_DIR)

  if(LPSolve_INCLUDE_DIR AND LPSolve_LIBRARIES)
    set(LPSolve_FOUND 1 CACHE BOOL "Found LPSolve library")
    message(STATUS "Found LPSolve Library: ${LPSolve_LIBRARIES}")
  else(LPSolve_INCLUDE_DIR AND LPSolve_LIBRARIES)
    set(LPSolve_FOUND 0 CACHE BOOL "Not found LPSolve library")
    message(STATUS "LPSolve Library: not found")
  endif(LPSolve_INCLUDE_DIR AND LPSolve_LIBRARIES)
endif(NOT LPSolve_INCLUDE_DIR)

mark_as_advanced(LPSolve_INCLUDE_DIR LPSolve_LIBRARIES LPSolve_DLL_DIR LPSolve_DLLS)
