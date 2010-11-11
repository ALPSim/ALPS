#!/bin/sh
#  Copyright Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_

cd $1
for f in *py
do
  vispython $f
done
