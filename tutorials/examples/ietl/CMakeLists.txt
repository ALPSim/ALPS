#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

#
# examples using the IETL
#

set(EXAMPLES
    arnoldi1_real
    arnoldi1_complex
    lanczos1
    lanczos2
    lanczos3
    lanczos4
    lanczos5
    power1
    power2
    jacobidavidson1
    jacobidavidson2
    lanczos_checkpoint
    lanczos_resume
    lanczos_checkpoint_complex
    lanczos_resume_complex
    jacobi) # inverse1 inverse2 currently fail jacobi also
foreach(program ${EXAMPLES})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} PRIVATE ALPS::alps)
endforeach()
