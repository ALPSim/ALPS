# SPDX-License-Identifier: MIT
function(alps_install_compiler_runtime destination)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
  include(InstallRequiredSystemLibraries)
  # MSVC's ARM64 redist also contains an ARM64EC-only exception runtime.
  # CMake includes it by filename, although native ARM64 cannot load it.
  if(CMAKE_CXX_COMPILER_ARCHITECTURE_ID STREQUAL "ARM64")
    list(FILTER CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS EXCLUDE REGEX "/vcruntime140_1d?\\.dll$")
  endif()
  install(
    PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
    DESTINATION "${destination}"
    COMPONENT runtime)
  install(
    DIRECTORY ${CMAKE_INSTALL_SYSTEM_RUNTIME_DIRECTORIES}
    DESTINATION "${destination}"
    COMPONENT runtime)
endfunction()
