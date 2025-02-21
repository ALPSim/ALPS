#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

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

