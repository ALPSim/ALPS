#!/bin/sh

branch=bindings
taghead=BINDINGS
checkfile=lapack/syev.hpp

dir=src/boost/numeric/bindings
files='blas lapack traits'
import_option="-m '' -ko"

srcdir=$1
tmpdir=import.$$

this=`basename $0`
if test -f "$this"; then :; else
  echo "ERROR: this script should be run in the config directory" && exit -1
fi
if test -z "$srcdir"; then
  echo "ERROR: $0 src_dir" && exit -1
fi
if test -f "../$dir/$checkfile"; then :; else
  echo "ERROR: ../$dir/$checkfile not found" && exit -1
fi
if test -f "$srcdir/$checkfile"; then :; else
  echo "ERROR: import source files ($checkfile) not found" && exit -1
fi

ROOT=`cat ../$dir/CVS/Root`
REPOSITORY=`cat ../$dir/CVS/Repository`
echo "CVS root = $ROOT"
echo "CVS repository = $REPOSITORY"

OLD_TAG=`(cd ../$dir && cvs status -v $checkfile | grep $taghead | head -1 | awk '{print $1}')`
echo -n "Tag for previously imported source files [$OLD_TAG] ? : "
read answer &> /dev/null
if test -n "$answer"; then
  OLD_TAG="$answer"
fi
if test -z "$OLD_TAG"; then
  echo "ERROR: empty tag." && exit -1
fi

NEW_TAG=${taghead}_`date +'%Y%m%d'`
echo -n "Tag for newly imported source files [$NEW_TAG] ? : "
read answer &> /dev/null
if test -n "$answer"; then
  NEW_TAG="$answer"
fi

if test "$NEW_TAG" = "$OLD_TAG"; then
  echo "ERROR: tag $NEW_TAG already exists" && exit -1
fi

mkdir $tmpdir

echo "Copying files into $tmpdir..."
for i in $files; do
  cp -rp $srcdir/$i $tmpdir/
done

set -x

(cd $tmpdir && cvs -d ${ROOT} import ${import_option} ${REPOSITORY} ${branch} ${NEW_TAG})

(cd ../$dir && cvs up -j $OLD_TAG -j $NEW_TAG)

set +x

rm -rf $tmpdir
