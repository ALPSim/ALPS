#  Copyright Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


install(FILES FindBoostSrc.cmake FindDocbookDtd.cmake FindDocbookXsl.cmake
              FindFop.cmake FindLPSolve.cmake FindVistrails.cmake
              FindXsltproc.cmake
        DESTINATION share/cmake COMPONENT build)
