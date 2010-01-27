#!/bin/sh
#  Copyright Synge Todo 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

taghead=BINDINGS
checkfile=lapack/syev.hpp
dir=src/boost/numeric/bindings
files='blas lapack traits'

this=`basename $0`
if test -f "$this"; then :; else
  echo "ERROR: this script should be run in the config directory" && exit -1
fi
if test -f "../$dir/$checkfile"; then :; else
  echo "ERROR: ../$dir/$checkfile not found" && exit -1
fi

TAG=`(cd ../$dir && cvs status -v $checkfile | awk 'NF==3' | grep $taghead | head -1 | awk '{print $1}')`

set -x

(cd ../$dir && cvs diff -r ${TAG} $files)
