#!/bin/sh
#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

mkdir -p deplibs/lib
mkdir -p deplibs/bin
cp /opt/alps/lib/libhdf5*dylib deplibs/lib
cp /opt/alps/lib/libhdf5*a deplibs/lib
cp /opt/local/lib/libsz*dylib deplibs/lib
cp /opt/local/lib/libz*dylib deplibs/lib
cp /opt/alps/lib/libmpi* deplibs/lib
cp /opt/alps/lib/libvt* deplibs/lib
cp /opt/alps/lib/libotf* deplibs/lib
cp /opt/alps/lib/libmca_common_sm* deplibs/lib
cp /opt/alps/lib/libopen-* deplibs/lib
cp /opt/alps/lib/libgfortran* deplibs/lib
cp /opt/alps/lib/libgcc* deplibs/lib
cp /opt/alps/lib/libstdc+* deplibs/lib
cp /opt/alps/lib/libquadmath* deplibs/lib
cp /opt/alps/lib/libgomp* deplibs/lib
cp -r /opt/alps/lib/openmpi deplibs/lib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.dylib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.2.dylib
install_name_tool -id /opt/alps/lib/libsz.2.dylib deplibs/lib/libsz.2.0.0.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.1.dylib
install_name_tool -id /opt/alps/lib/libz.1.dylib deplibs/lib/libz.1.2.*.dylib
install_name_tool -id /opt/alps/lib/libzmq.5.dylib deplibs/lib/libzmq.5.dylib
#install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5.dylib
#install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5.*.dylib
#install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5_hl.dylib
#install_name_tool -change /opt/local/lib/libsz.2.dylib /opt/alps/lib/libsz.2.dylib deplibs/lib/libhdf5_hl.6.dylib
#install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5.dylib
#install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5.6.dylib
#install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5_hl.dylib
#install_name_tool -change /opt/local/lib/libz.1.dylib /opt/alps/lib/libz.1.dylib deplibs/lib/libhdf5_hl.6.dylib

mkdir -p deplibs/include
cp /opt/alps/include/H5* deplibs/include
cp /opt/alps/include/hdf* deplibs/include
cp /opt/local/include/zconf.h deplibs/include
cp /opt/local/include/zlib.h deplibs/include
cp /opt/local/include/szlib.h deplibs/include
cp /opt/local/include/szip_adpt.h deplibs/include
cp /opt/local/include/ricehdf.h deplibs/include
cp /opt/alps/include/mpi.h deplibs/include
cp -r /opt/alps/include/openmpi deplibs/include
cp -r /opt/alps/include/vampirtrace deplibs/include

mkdir -p deplibs/bin
cp /opt/alps/bin/mpi* deplibs/bin
cp /opt/alps/bin/ompi* deplibs/bin
cp /opt/alps/bin/orte* deplibs/bin
cp /opt/alps/bin/opal_wrapper deplibs/bin

mkdir -p deplibs/etc
cp /opt/alps/etc/openmpi* deplibs/etc

mkdir -p deplibs/share
cp -r /opt/alps/share/openmpi deplibs/share
cp -r /opt/alps/share/man deplibs/share

sudo cp -r deplibs/lib/* /opt/alps/lib
sudo cp -r deplibs/include/* /opt/alps/include
sudo cp deplibs/bin/* /opt/alps/bin
sudo cp deplibs/etc/* /opt/alps/etc
sudo cp -r deplibs/share/* /opt/alps/share
