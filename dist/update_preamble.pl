#!/usr/local/bin/perl

# Script for updating preamble of *.h and *.C
# by Synge Todo <wistaria@comp-phys.org>

$basedir = $0;
$basedir =~ s/[a-zA-Z\_\.]+$//g;
$skel = join('', $basedir, "preamble.in");

foreach $file (@ARGV) {
    if (-f $file) {
	$file_new = "$file.new";

	# scan
	$year0 = "";
	$year1 = "";
	$id = "";
	@authors = ();
	@emails = ();
	$finish_preamble = 0;
	$skip = 0;
	open(ORIG, "< $file");
	open(NEW, "> $file.new");
	foreach $line (<ORIG>) {
	    $line =~ s/\t/        /g;
	    chomp($line);
	    if ($finish_preamble == 0) {
		if ($line =~ /^\*\s+Copyright.+([0-9]{4})-([0-9]{4})\s+by\s+(\S\C+)\s+\<([a-zA-Z0-9\.\-_]+@[a-zA-Z0-9\.\-_]+)\>/) {
		    $year0 = $1;
		    $year1 = $2;
		    @authors[$#authors+1] = $3;
		    @emails[$#emails+1] = $4;
		} elsif ($line =~ /^\*\s+Copyright.+([0-9]{4})\s+by\s+(\S\C+)\s+\<([a-zA-Z0-9\.\-_]+@[a-zA-Z0-9\.\-_]+)\>/) {
		    $year0 = $1;
		    @authors[$#authors+1] = $2;
		    @emails[$#emails+1] = $3;
		} elsif ($line =~ /^\*\s+(\S\C+)\s+\<([a-zA-Z0-9\.\-_]+@[a-zA-Z0-9\.\-_]+)\>/) {
		    @authors[$#authors + 1] = $1;
		    @emails[$#emails + 1] = $2;
		} elsif ($line =~ /(\$Id\:\C+\$)/) {
		    $id = $1;
		} elsif ($line =~ /(^[\s\*\/])|(^$)/) {
		    ## nothing to do
		} else {
		    $finish_preamble = 1;
		}
		if ($finish_preamble == 1) {

		    ## check
		    if ($year0 eq "" || @authors[0] eq "") {
			$skip = 1;
		    } else {
			if ($year1 eq $year0) { $year1 = ""; }
			if ($id eq "") { $id = "\$Id$"; }
			
			## print out preamble
			open(SKEL, "< $skel") || die "Couldn't open $skel";
			foreach $sk (<SKEL>) {
			    chomp($sk);
			    if ($sk =~ /\@COPYRIGHT\@/) {
				if ($year1) {
				    print NEW "* Copyright (C) $year0-$year1 by @authors[0] <@emails[0]>";
				} else {
				    print NEW "* Copyright (C) $year0 by @authors[0] <@emails[0]>";
				}
				if ($#authors > 0) { print NEW ","; }
				print NEW "\n";
				for ($i = 1; $i <= $#authors; $i++) {
				    if ($year1) {
					print NEW "*                            @authors[$i] <@emails[$i]>";
				    } else {
					print NEW "*                       @authors[$i] <@emails[$i]>";
				    }
				    if ($i < $#authors) { print NEW ","; }
				    print NEW "\n";
				}
			    } elsif ($sk =~ /\@ID\@/) {
				print NEW "// $id\n";
			    } else {
				print NEW "$sk\n";
			    }
			}
			print NEW "\n";
		    }
		}
	    }

	    if ($skip == 0 && $finish_preamble == 1) {
		print NEW "$line\n";
	    }
	}
	close(ORIG);
	close(NEW);

	if ($skip ==0) {
	    system("diff $file $file_new > /dev/null");
	    if ($? == 256) {
		unlink $file;
		rename $file_new, $file;
		print "$file is updated.\n";
	    } else {
		unlink $file_new;
	    }
	} else {
	    print "$file does not obey ALPS standard.  Skipped.\n";
	    unlink $file_new;
	}
    }
}
