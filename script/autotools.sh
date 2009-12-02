#!/bin/sh

AC_VERSION=2.65
AM_VERSION=1.11
LT_VERSION=2.2.6b

PREFIX="$1"
BUILD_DIR="$2"
SRC_DIR="$3"

if test -z "$BUILD_DIR"; then
  echo "$0 prefix build_dir [src_dir]"
  exit 127
fi

# for Mac OS X (which has apple version of libtool)
if test `libtool --version > /dev/null 2>&1; echo $?` != 0; then
  LT_OPT="--program-prefix=g"
fi

LOG="$0.log.$$"
echo "executing $0 $*" | tee "$LOG"

# autoconf

echo "building autoconf..." | tee -a "$LOG"

URL="http://ftp.gnu.org/gnu/autoconf/autoconf-$AC_VERSION.tar.gz"
SRC="$SRC_DIR/autoconf-$AC_VERSION.tar.gz"

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/autoconf-$AC_VERSION"
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

( \
echo "configuring..." && \
(cd "$BUILD_DIR/autoconf-$AC_VERSION" && ./configure --prefix="$PREFIX") && \
echo "building..." && \
(cd "$BUILD_DIR/autoconf-$AC_VERSION" && make) && \
echo "installing..." && \
(cd "$BUILD_DIR/autoconf-$AC_VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/autoconf-$AC_VERSION" \
) 2>&1 | tee -a "$LOG"

# automake

echo "building automake..." | tee -a "$LOG"

URL="http://ftp.gnu.org/gnu/automake/automake-$AM_VERSION.tar.gz"
SRC="$SRC_DIR/automake-$AM_VERSION.tar.gz"

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/automake-$AM_VERSION"
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

( \
echo "configuring..." && \
(cd "$BUILD_DIR/automake-$AM_VERSION" && ./configure --prefix="$PREFIX") && \
echo "building..." && \
(cd "$BUILD_DIR/automake-$AM_VERSION" && make) && \
echo "installing..." && \
(cd "$BUILD_DIR/automake-$AM_VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/automake-$AM_VERSION" \
) 2>&1 | tee -a "$LOG"

# libtool

echo "building libtool..." | tee -a "$LOG"

URL="http://ftp.gnu.org/gnu/libtool/libtool-$LT_VERSION.tar.gz"
SRC="$SRC_DIR/libtool-$LT_VERSION.tar.gz"

echo "cleaning up..." | tee -a "$LOG"
if test -d "$BUILD_DIR"; then
  rm -rf "$BUILD_DIR/libtool-$LT_VERSION"
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

( \
echo "configuring..." && \
(cd "$BUILD_DIR/libtool-$LT_VERSION" && ./configure --prefix="$PREFIX" "$LT_OPT") && \
echo "building..." && \
(cd "$BUILD_DIR/libtool-$LT_VERSION" && make) && \
echo "installing..." && \
(cd "$BUILD_DIR/libtool-$LT_VERSION" && make install) && \
echo "cleaning up..." && \
rm -rf "$BUILD_DIR/libtool-$LT_VERSION" \
) 2>&1 | tee -a "$LOG"

echo "done" | tee -a "$LOG"
echo "log file = $LOG" | tee -a "$LOG"
