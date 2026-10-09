#!/bin/sh
#  Copyright Synge Todo and Matthias Troyer 2009 - 2011.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

VERSION=4.9.3
GMPVERSION=5.1.3
MPFRVERSION=3.1.3
MPCVERSION=1.0.3

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

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/gcc-$VERSION"
  rm -rf "$BUILD_DIR/gmp-$GMPVERSION"
  rm -rf "$BUILD_DIR//mpfr-$MPFRVERSION"
  rm -rf "$BUILD_DIR/mpc-$MPCVERSION"
  rm -rf "$BUILD_DIR/aux"
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
(cd "$BUILD_DIR" && curl "$GMPURL" } bzcat | tar xf -  && curl "$MPFRURL" | tar zxf - && curl "$MPCURL" | tar zxf -) 2>&1 | tee -a "$LOG"

echo "building gmp ... " && \
(cd "$BUILD_DIR/gmp-$GMPVERSION" && ./configure --prefix=$PREFIX && make install) 2>&1 | tee -a "$LOG"

echo "building mpfr ... " && \
(cd "$BUILD_DIR/mpfr-$MPFRVERSION" && ./configure --prefix=$PREFIX --with-gmp=$PREFIX && make install) 2>&1 | tee -a "$LOG"

echo "building mpc ... " && \ 
(cd "$BUILD_DIR/mpc-$MPCVERSION" && ./configure --prefix=$PREFIX --with-gmp=$PREFIX --with-mpfr=$PREFIX && make install) 2>&1 | tee -a "$LOG"
( \
echo "configuring..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && ./configure --enable-languages=fortran --prefix="$PREFIX" --with-gmp=$PREFIX --with-mpfr=$PREFIX --with-mpc=$PREFIX ) && \
echo "building..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && make bootstrap) && \
echo "installing..." && \
(cd "$BUILD_DIR/gcc-$VERSION" && make install) && \

echo "done" \
) 2>&1 | tee -a "$LOG"

echo "log file = $LOG" | tee -a "$LOG"
