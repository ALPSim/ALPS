#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

######################################################################
# applications
######################################################################

add_definitions(${LAPACK_DEFINITIONS})
add_subdirectory(mc)
add_subdirectory(qmc)
add_subdirectory(diag)
add_subdirectory(dmrg)
add_subdirectory(dmft)
  

######################################################################
# GUI's
######################################################################

# add_subdirectory(gui/simulator)

