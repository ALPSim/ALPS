#! /usr/bin/env python

import sys

ns = sys.argv[1]
nedges = str(int(ns) - 1)

print '<LATTICES>'
print ' <GRAPH name = "open chain lattice with special edges" dimension="1" vertices="' + ns + '" edges="' + nedges + '">'
print '  <VERTEX id="1" type="0"><COORDINATE>0</COORDINATE></VERTEX>'
for i in range(2,int(ns)):
    print '  <VERTEX id="' + str(i) + '" type="1"><COORDINATE>' + str(i) + '</COORDINATE></VERTEX>'
print ' <VERTEX id="' + ns + '" type="0"><COORDINATE>' + ns + '</COORDINATE></VERTEX>'
for i in range(1,int(ns)):
    nn = str(i + 1)
    print '  <EDGE source="' + str(i) + '" target="' + nn + '" id="' + str(i) + '" type="0" vector="1"/>'

print ' </GRAPH>'
print '</LATTICES>'
