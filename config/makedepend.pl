#!/usr/local/bin/perl
#
# PALM++/palm library
#
# makedepend.pl : supplemental perl script processing output of makedepend
#
# $Id$
#
# Copyright (C) 2002 by Synge Todo <wistaria@comp-phys.org>,
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# usage:
#   rm -f Makefile.dep && touch Makefile.dep
#   makedepend -f Makefile.dep -- -Y -I../../src -- *.C > /dev/null 2>&1
#   mv -f Makefile.dep makedepend.tmp
#   cat makedepend.tmp | perl makedepend.pl > Makefile.dep
#   rm -f Makefile.dep.bak makedepend.tmp

$SRCDIR = '$(srcdir)';
$TOP_SRCDIR = '$(top_srcdir)';
$TOP_BUILDDIR = '$(top_builddir)';

while (<>) {
    chomp;
    @elems = split(' ');
    if (@elems[0] ne '#' && @elems[0] ne '') {
	$base = shift @elems;
	$base =~ s/\.o\://;
	foreach $header (@elems) {
	    # remove preceding '../'
	    $header =~ s/(\.\.\/)//g;
	    if ($header =~ /^src/ || $header =~ /^include/) {
		if ($header =~ /config\.h$/) {
		    $header = "$TOP_BUILDDIR/$header";
		} else {
		    $header = "$TOP_SRCDIR/$header";
		}
	    } else {
		$header = "$SRCDIR/$header";
	    }
	    $headers{$base} = "$headers{$base} $header";
	}
    }
}

foreach $base (keys(%headers)) {
    print "${base}.o ${base}_mpi.o ${base}_pvm.o ${base} ${base}_mpi ${base}_pvm:$headers{$base}\n";
}
