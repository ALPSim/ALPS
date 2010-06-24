#!/usr/bin/perl
#
# Generate a lattice file for the spin S=1 chain with especial S=1/2 edges 
#

$ns = $ARGV[0];
$nedges = $ns - 1;
print "<LATTICES>\n";
print " <GRAPH name = \"open chain lattice with special edges\" dimension=\"1\" vertices=\" $ns \" edges=\"$nedges\">\n";
print "  <VERTEX id=\"1\" type=\"0\"><COORDINATE>0</COORDINATE></VERTEX>\n";
for($i = 2; $i < $ns; $i = $i + 1)
{
  print "  <VERTEX id=\"$i\" type=\"1\"><COORDINATE>$i</COORDINATE></VERTEX>\n";
}
print " <VERTEX id=\"$ns\" type=\"0\"><COORDINATE>$ns</COORDINATE></VERTEX>\n";
for($i = 1; $i < $ns; $i = $i + 1)
{
 $nn = $i + 1;
 print "  <EDGE source=\"$i\" target=\"$nn\" id=\"$i\" type=\"0\" vector=\"1\"/>\n";
}

print " </GRAPH>\n";
print "</LATTICES>\n";
