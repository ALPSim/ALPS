# - Find SQLite
# Find SQLite headers and libraries.
#
#  SQLite_FOUND        - True if SQLite found.
#  SQLite_INCLUDE_DIR  - where to find sqlite3.h
#  SQLite_LIBRARIES    - list of libraries when using SQLite.
#  SQLite_DLLS         - list of DLLs when using SQLite.

if(NOT SQLite_INCLUDE_DIR)
  set(__TRIAL_PATHS "$ENV{SQLITE_HOME}/include" /usr/include
      /usr/local/include /opt/include
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/include"
     )
  find_path(SQLite_INCLUDE_DIR sqlite3.h ${__TRIAL_PATHS})
  set(__TRIAL_PATHS "$ENV{SQLITE_HOME}/lib" /usr/lib
      /usr/local/lib /opt/lib
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/lib"
     )
  find_library(SQLite_LIBRARIES sqlite3 ${__TRIAL_PATHS})
  set(__TRIAL_PATHS "$ENV{SQLITE_HOME}/bin"
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/opt/bin"
     )
  find_path(SQLite_DLL_DIR sqlite3.dll ${__TRIAL_PATHS})
  if(SQLite_DLL_DIR)
    set(SQLite_DLLS "${SQLite_DLL_DIR}/sqlite3.dll" CACHE PATH "SQLite DLL")
  endif(SQLite_DLL_DIR)

  if(SQLite_INCLUDE_DIR AND SQLite_LIBRARIES)
    set(SQLite_FOUND 1 CACHE BOOL "Found SQLite library")
    message(STATUS "Found SQLite Library: ${SQLite_LIBRARIES}")
  else(SQLite_INCLUDE_DIR AND SQLite_LIBRARIES)
    set(SQLite_FOUND 0 CACHE BOOL "Not found SQLite library")
    message(STATUS "SQLite Library: not found")
  endif(SQLite_INCLUDE_DIR AND SQLite_LIBRARIES)
endif(NOT SQLite_INCLUDE_DIR)

mark_as_advanced(SQLite_INCLUDE_DIR SQLite_LIBRARIES SQLite_DLL_DIR SQLite_DLLS)
