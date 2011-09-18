#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# examples using the IETL
# 

if (ALPS_BUILD_EXAMPLES)
  if(LAPACK_FOUND)
    add_definitions(${LAPACK_DEFINITIONS})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
    set (EXAMPLES arnoldi1_real arnoldi1_complex lanczos1 lanczos2 lanczos3 lanczos4 lanczos5 power1 power2 jacobidavidson1 jacobidavidson2 lanczos_checkpoint lanczos_resume lanczos_checkpoint_complex lanczos_resume_complex jacobi) # inverse1 inverse2 currently fail jacobi also
    foreach(program ${EXAMPLES})
      add_executable(${program} ${program}.cpp)
      target_link_libraries(${program} alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
    endforeach(program)
  endif(LAPACK_FOUND)
endif (ALPS_BUILD_EXAMPLES)

