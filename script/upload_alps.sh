#!/bin/sh
#  Copyright Matthias Troyer 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


scp alps-vistrails*.*[gez] root@alps.ethz.ch:/var/www/vistrails
scp alps-download*.*[gez] root@alps.ethz.ch:/var/www/software/releases
scp alps-src*.*z root@alps.ethz.ch:/var/www/software/releases
scp alps-2*.*[gez] root@alps.ethz.ch:/var/www/software/releases
scp CPackUploads/*pkg CPackUploads/*zip root@alps.ethz.ch:/var/www/software/releases/packages
