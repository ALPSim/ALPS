#!/bin/sh
#  Copyright Synge Todo and Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

echo "ALPS Libraries:"
echo C++:; find legacy lib src tool vistrails -name "*.h" -or -name "*.hpp" -or -name "*.C" -or -name "*.cpp" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo py:; find legacy lib src tool vistrails -name "*.py" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo f90:; find legacy lib src tool vistrails -name "*.f90" | grep -v svn | grep -v boost | xargs wc -l | tail -1

echo "ALPS Applications:"
echo C++:; find applications -name "*.h" -or -name "*.hpp" -or -name "*.C" -or -name "*.cpp" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo py:; find applications -name "*.py" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo f90:; find applications -name "*.f90" | grep -v svn | grep -v boost | xargs wc -l | tail -1

echo "Examples, tutorials, tests:"
echo C++:; find example tutorials test -name "*.h" -or -name "*.hpp" -or -name "*.C" -or -name "*.cpp" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo py:; find example tutorials test -name "*.py" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo f90:; find example tutorials test -name "*.f90" | grep -v svn | grep -v boost | xargs wc -l | tail -1
