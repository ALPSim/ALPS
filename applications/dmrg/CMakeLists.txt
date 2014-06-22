#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

option(ALPS_BUILD_MPS "Build and install MPS application" ON)

if(ALPS_BUILD_MPS)
  add_subdirectory(mps)
endif(ALPS_BUILD_MPS)
add_subdirectory(dmrg)
add_subdirectory(tebd)
