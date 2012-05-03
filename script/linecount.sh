#!/bin/sh
#  Copyright Synge Todo and Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
echo h:; find . -name "*.h" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo hpp:; find . -name "*.hpp" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo C:; find . -name "*.C" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo cpp:; find . -name "*.cpp" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo py:; find . -name "*.py" | grep -v svn | grep -v boost | xargs wc -l | tail -1
echo f90:; find . -name "*.f90" | grep -v svn | grep -v boost | xargs wc -l | tail -1
