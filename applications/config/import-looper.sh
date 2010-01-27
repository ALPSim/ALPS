#!/bin/sh
#  Copyright Synge Todo 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# for echo without trailing newline character
case `echo "testing\c"; echo 1,2,3`,`echo -n testing; echo 1,2,3` in
  *c*,-n*) ECHO_N= ECHO_C='
' ;;
  *c*,*  ) ECHO_N=-n ECHO_C= ;;
  *)       ECHO_N= ECHO_C='\c' ;;
esac

checkfile=loop.C
dir=qmc/looper
files_remove="CITATIONS.txt LICENSE.txt Makefile.am Makefile.in README acinclude.m4 aclocal.m4 bootstrap check config configure configure.ac diag.C doc extras ising.C loop_ns* qwl_evaluate.C setup.sh sse_qwl.C standalone test v3.1 looper/correlation_length.h looper/find_bridge.h looper/flatten_matrix.h looper/gap.h looper/generate_seed.h looper/histogram.h looper/localsus.h looper/parallel.h looper/poisson_distribution.h looper/sop.h looper/top.h looper/transmag.h"
convert_preamble=yes
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

REPOSITORY=`(cd $srcdir && svn info) | grep URL | awk '{print $2}'`

echo "SVN repository = $REPOSITORY"

rm -rf $tmpdir

echo "Checking out files into $tmpdir..."
svn co $REPOSITORY $tmpdir
for i in $files_remove; do
  rm -rf $tmpdir/$i
done

REVISION=`svn info $tmpdir | grep 'Last Changed Rev:' | awk '{print $4}'`

if test "$convert_preamble" = yes; then
  echo "Converting preambles..."
  find $tmpdir -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs ./update_preamble
fi

find $tmpdir -name '.svn' | xargs rm -rf

echo "SVN version = $REVISION"

echo "repository: $REPOSITORY" > ../$dir/IMPORT
echo "revirsion: $REVISION" >> ../$dir/IMPORT

cp -r $tmpdir/* ../$dir

rm -rf $tmpdir
