#!/bin/sh
#  Copyright Synge Todo 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if test -z "$1"; then
  echo "Usage: $0 dir"
  exit -1
fi

DIR=$1
if test -d "$DIR"; then 
  echo "Downloading files to $DIR ..."
else
  echo "Directory \"$DIR\" does not exists.  Stop."
  exit -1
fi

# User configuration
DOCBOOK_XSL_VERSION=1.69.1
DOCBOOK_XSL_TARBALL=docbook-xsl-$DOCBOOK_XSL_VERSION.tar.gz
DOCBOOK_XSL_URL=http://puzzle.dl.sourceforge.net/sourceforge/docbook/$DOCBOOK_XSL_TARBALL
DOCBOOK_XSL_DIR=$DIR/docbook-xsl-$DOCBOOK_XSL_VERSION

DOCBOOK_DTD_VERSION=4.4
DOCBOOK_DTD_ZIP=docbook-xml-$DOCBOOK_DTD_VERSION.zip
DOCBOOK_DTD_URL=http://www.oasis-open.org/docbook/xml/$DOCBOOK_DTD_VERSION/$DOCBOOK_DTD_ZIP
DOCBOOK_DTD_DIR=$DIR/docbook-dtd-$DOCBOOK_DTD_VERSION

FOP_VERSION=0.20.5
FOP_TARBALL=fop-$FOP_VERSION-bin.tar.gz
FOP_URL=http://mirrors.ibiblio.org/pub/mirrors/apache/xml/fop/$FOP_TARBALL
FOP_DIR=$DIR/fop-$FOP_VERSION

# Get the DocBook XSLT Stylesheets
if test -d $DOCBOOK_XSL_DIR; then
  echo "DocBook XSLT Stylesheets already exists (version $DOCBOOK_XSL_VERSION).  Skip."
else
  echo "Downloading DocBook XSLT Stylesheets version $DOCBOOK_XSL_VERSION ..."
  curl $DOCBOOK_XSL_URL | (cd $DIR && tar zxf -)
fi

# Get the DocBook DTD
if test -d $DOCBOOK_DTD_DIR; then
  echo "DocBook XML DTD already exists (version $DOCBOOK_DTD_VERSION).  Skip."
else
  echo "Downloading DocBook XML DTD version $DOCBOOK_DTD_VERSION ..."
  (cd $DIR && curl -O $DOCBOOK_DTD_URL && unzip -q $DOCBOOK_DTD_ZIP -d $DOCBOOK_DTD_DIR && rm -f $DOCBOOK_DTD_ZIP)
fi

# Get Apache FOP
if test -d $FOP_DIR; then
  echo "Apache FOP already exists (version $FOP_VERSION).  Skip."
else
  echo "Downloading Apache FOP version $FOP_VERSION ..." 
  curl $FOP_URL | (cd $DIR && tar zxf -)
fi

echo "Done with everything."
