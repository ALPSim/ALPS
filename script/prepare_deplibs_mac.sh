#!/bin/sh
#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

mkdir -p deplibs/lib
mkdir -p deplibs/bin
cp /opt/local/lib/libhdf5*dylib deplibs/lib
cp /opt/local/lib/libsz*dylib deplibs/lib
cp /opt/local/lib/libz*dylib deplibs/lib
install_name_tool -id /opt/alps/lib/libhdf5.6.dylib deplibs/lib/libhdf5.dylib
install_name_tool -id /opt/alps/lib/libhdf5.6.dylib deplibs/lib/libhdf5.6.dylib
install_name_tool -id /opt/alps/lib/libhdf5_hl.6.dylib deplibs/lib/libhdf5_hl.dylib
install_name_tool -id /opt/alps/lib/libhdf5_hl.6.dylib deplibs/lib/libhdf5_hl.6.dylib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.dylib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.2.dylib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.2.0.0.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.1.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.1.2.*.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5.*.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5_hl.dylib
install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5_hl.6.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5.6.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5_hl.dylib
install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5_hl.6.dylib

mkdir -p deplibs/include
cp /opt/local/include/H5* deplibs/include
cp /opt/local/include/hdf* deplibs/include
cp /opt/local/include/zconf.h deplibs/include
cp /opt/local/include/zlib.h deplibs/include
cp /opt/local/include/szlib.h deplibs/include
cp /opt/local/include/szip_adpt.h deplibs/include
cp /opt/local/include/ricehdf.h deplibs/include

sudo cp deplibs/lib/* /opt/alps/lib
sudo cp -r deplibs/include/* /opt/alps/include
