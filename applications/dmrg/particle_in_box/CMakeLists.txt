#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

add_executable(simple_dmrg dmrgmain.C dmrg3.C)
target_link_libraries(simple_dmrg alps boost ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
install(TARGETS simple_dmrg RUNTIME DESTINATION bin COMPONENT applications)
