# CT-HYB-SEGMENT: segment-picture hybridization expansion CT-QMC solver,
# formerly ALPSCore/CT-HYB-SEGMENT (Hafermann, Werner, Gull, CPC 184, 1280 (2013)).
# SPDX-License-Identifier: MIT

# The solver reads and writes Green's functions through ALPS::gf, which needs Eigen.
if(NOT TARGET ALPS::gf)
  message(STATUS "ALPS::gf is not built (Eigen not found): ct-hyb-segment is not built")
  return()
endif()

set(_cthyb_segment_sources
    run.cpp
    hybsim.cpp
    hyblocal.cpp
    hybint.cpp
    hybfun.cpp
    hybretintfun.cpp
    hybmatrix.cpp
    hybconfig.cpp
    hybupdates.cpp
    hybevaluate.cpp
    hybmeasurements.cpp)

# NFFT3 evaluates the frequency measurements by a non-equispaced FFT; without it,
# hybmatrix_ft.cpp computes the same quantities by direct summation.
if(ALPS_ENABLE_NFFT)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(NFFT3 REQUIRED)
  list(APPEND _cthyb_segment_sources hybmatrix_nfft.cpp)
else()
  list(APPEND _cthyb_segment_sources hybmatrix_ft.cpp)
endif()

add_library(alps_cthyb_segment STATIC ${_cthyb_segment_sources})
add_library(ALPS::cthyb_segment ALIAS alps_cthyb_segment)
set_target_properties(alps_cthyb_segment PROPERTIES EXPORT_NAME cthyb_segment
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_link_libraries(alps_cthyb_segment PUBLIC ALPS::alps ALPS::gf)
target_include_directories(alps_cthyb_segment PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
# alps/gf/flagcheck.hpp asks consumers to choose its debug checks, as ALPS::gf does.
target_compile_definitions(alps_cthyb_segment PUBLIC
  "$<IF:$<CONFIG:Debug>,ALPS_GF_DEBUG,BOOST_DISABLE_ASSERTS>")
if(ALPS_ENABLE_NFFT)
  target_include_directories(alps_cthyb_segment PRIVATE ${NFFT3_INCLUDE_DIRS})
  target_link_libraries(alps_cthyb_segment PRIVATE ${NFFT3_LIBRARIES})
endif()

add_executable(ct-hyb-segment hybmain.cpp)
target_link_libraries(ct-hyb-segment PRIVATE ALPS::cthyb_segment)
alps_install_applications(ct-hyb-segment)
