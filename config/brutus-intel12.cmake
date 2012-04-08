#  Copyright Jan Gukelberger 2012.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# We assume the following modules are loaded:
#   module load hdf5 goto2 python/2.7.2 intel/12 open_mpi/1.5 cmake
# Do NOT load module mkl, it is already included with intel/12.
#
# icpc option "-fast" enables "-static" which does not work well with cmake's
# library search when there are dynamic libraries around. 
# Therefore, in CMAKE_CXX_FLAGS_RELEASE we manually set allother options
# included in "-fast". 

SET(CMAKE_CXX_COMPILER $ENV{CXX})
SET(CMAKE_C_COMPILER $ENV{CC})
SET(CMAKE_Fortran_COMPILER $ENV{FC})
SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost -mkl -DNDEBUG -DBOOST_DISABLE_ASSERTS" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -mkl" CACHE STRING "" FORCE)
SET(CMAKE_AR "/cluster/apps/intel/composerxe-2011.3.174/bin/intel64/xiar" CACHE STRING "" FORCE) 
SET(CMAKE_LINKER "/cluster/apps/intel/composerxe-2011.3.174/bin/intel64/xild" CACHE STRING "" FORCE)
SET(LAPACK_FOUND 1) # with -mkl we have lapack no matter what cmake thinks
