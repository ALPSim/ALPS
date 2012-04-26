#!/bin/sh
#  Copyright Synge Todo and Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

VERSION=4.4.7
GMPVERSION=5.0.4
MPFRVERSION=3.1.0
MPCVERSION=0.8.2

PREFIX="$1"
BUILD_DIR="$2"
SRC_DIR="$3"

if test -z "$BUILD_DIR"; then
  echo "$0 prefix build_dir [src_dir]"
  exit 127
fi

LOG="$0.log.$$"
echo "executing $0 $*" | tee "$LOG"

URL="ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/gcc-$VERSION/gcc-$VERSION.tar.gz"
GMPURL=ftp://ftp.gnu.org/gnu/gmp/gmp-$GMPVERSION.tar.bz2
MPFRURL=http://www.mpfr.org/mpfr-current/mpfr-$MPFRVERSION.tar.gz
MPCURL=http://www.multiprecision.org/mpc/download/mpc-$MPCVERSION.tar.gz
SRC="$SRC_DIR/gcc-$VERSION.tar.gz"
BDIR=$BUILD_DIR

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/gcc-$VERSION"
  rm -rf "$BUILD_DIR/gmp-$GMPVERSION"
  rm -rf "$BUILD_DIR//mpfr-$MPFRVERSION"
  rm -rf "$BUILD_DIR/mpc-$MPCVERSION"
  rm -rf "$BDIR/aux"
else
  mkdir -p "$BUILD_DIR"
fi

echo "retrieving source files..." | tee -a "$LOG"
if test -n "$SRC_DIR" && test -f "$SRC"; then
  (cd "$BUILD_DIR" && tar zxf $SRC) 2>&1 | tee -a "$LOG"
else
  CURL=`which curl`
  if test -z "$CURL"; then
    echo "curl utility not found" | tee -a "$LOG"
    exit 127
  fi
  (cd "$BUILD_DIR" && "$CURL" "$URL" | tar zxf -) 2>&1 | tee -a "$LOG"
fi
echo "downloading more source files" && \
(cd "$BUILD_DIR" && curl "$GMPURL" | bzcat | tar xf -  && curl "$MPFRURL" | tar zxf - && curl "$MPCURL" | tar zxf -) 2>&1 | tee -a "$LOG"

echo "building gmp ... " && \
(cd "$BUILD_DIR/gmp-$GMPVERSION" && ./configure --prefix=$BDIR/aux && make install) 2>&1 | tee -a "$LOG"

echo "building mpfr ... " && \
(cd "$BUILD_DIR/mpfr-$MPFRVERSION" && ./configure --prefix=$BDIR/aux --with-gmp=$BDIR/aux && make install) 2>&1 | tee -a "$LOG"

echo "building mpc ... " && \ 
(cd "$BUILD_DIR/mpc-$MPCVERSION" && ./configure --prefix=$BDIR/aux --with-gmp=$BDIR/aux --with-mpfr=$BDIR/aux && make install) 2>&1 | tee -a "$LOG"
( \
echo "configuring..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && ./configure --prefix="$PREFIX" --with-gmp=$BDIR/aux --with-mpfr=$BDIR/aux --with-mpc=$BDIR/aux ) && \
echo "building..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && make bootstrap) && \
echo "installing..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/gcc-$VERSION" && \
rm -rf "$BUILD_DIR/gmp-$GMPVERSION" && \
rm -rf "$BUILD_DIR//mpfr-$MPFRVERSION" && \
rm -rf "$BUILD_DIR/mpc-$MPCVERSION" && \
rm -rf "$BDIR/aux" && \

echo "done" \
) 2>&1 | tee -a "$LOG"

echo "log file = $LOG" | tee -a "$LOG"
