#!/usr/local/bin/perl

# makedepend.pl : supplemental perl script processing output of makedepend
#
# usage:
#   rm -f Makefile.dep && touch Makefile.dep
#   makedepend -f Makefile.dep -- -Y -I../../src -- `find * -name '*.C' -o -name '*.cpp'` > /dev/null 2>&1
#   mv -f Makefile.dep makedepend.tmp
#   cat makedepend.tmp | perl makedepend.pl | sort > Makefile.dep
#   rm -f Makefile.dep.bak makedepend.tmp

# written by Synge Todo <wistaria@comp-phys.org>
# $Id$

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
	    if ($header =~ /config\.h$/) {
		$header = "$TOP_BUILDDIR/$header";
	    } else {
		if ($header =~ /^src/ || $header =~ /^include/) {
		    $header = "$TOP_SRCDIR/$header";
		} else {
		    if ($header =~ /\//) {
			$header = "$TOP_SRCDIR/$header";
		    } else {
			$header = "$SRCDIR/$header";
		    }
		}
	    }
	    # remove './'
	    $header =~ s/(\/\.\/)/\//g;
	    $headers{$base} = "$headers{$base} $header";
	}
    }
}

foreach $base (keys(%headers)) {
    print "${base}.o ${base}_mpi.o ${base}_pvm.o ${base} ${base}_mpi ${base}_pvm:$headers{$base}\n";
}
