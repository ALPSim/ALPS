# $Id$

MAINTAINERCLEANFILES	= Makefile.in
dist_noinst_DATA = energy.xml lattices.xml susceptibility.xml \
	result.heisenberg_chain.parm.txt.output \
	result.heisenberg_ladder.parm.txt.output \
	result.ising_cluster.parm.txt.output \
	result.ising_local.parm.txt.output

.PHONY: test tests check

test tests check:
	-/bin/sh -e $(top_builddir)/test/mc/spins/check.sh ising_local.parm energy.xml
	-/bin/sh -e $(top_builddir)/test/mc/spins/check.sh ising_cluster.parm energy.xml
	-/bin/sh -e $(top_builddir)/test/mc/spins/check.sh heisenberg_chain.parm susceptibility.xml
	-/bin/sh -e $(top_builddir)/test/mc/spins/check.sh heisenberg_ladder.parm susceptibility.xml

clean-local:
	rm -rf *.out.xml
	rm -rf *.in.xml
	rm -rf *.out.run1
	rm -rf result.*.txt
