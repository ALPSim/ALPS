SUBDIRS =
if EXIST_DIAG
SUBDIRS += diag
endif
if EXIST_DMRG
SUBDIRS += dmrg
endif
if EXIST_DMFT
SUBDIRS += dmft
endif
if EXIST_MC
SUBDIRS += mc
endif
if EXIST_QMC
SUBDIRS += qmc
endif
if EXIST_TEST
SUBDIRS += test
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

REMOVE_FOR_RELEASE = dmrg/interacting qmc/sse-loop test/sse-loop  dmft

distdir = alps-applications-$(VERSION)

dist_noinst_SCRIPTS = config/update_preamble
dist_noinst_DATA = CITATIONS.txt ChangeLog LICENSE.txt config/preamble.in doc/index.html


dist_noinst_DATA += \
tutorial/lattices.xml tutorial/models.xml \
tutorial/classical1/parm1a tutorial/classical1/parm1b \
tutorial/classical2/parm2a tutorial/classical2/parm2b tutorial/classical2/plot.xml \
tutorial/quantum1/parm3a tutorial/quantum1/parm3b tutorial/quantum1/plot.xml \
tutorial/quantum2/parm4a tutorial/quantum2/parm4b tutorial/quantum2/parm4c tutorial/quantum2/plot2.xml \
tutorial/quantum3/parm5a tutorial/quantum3/parm5b \
tutorial/quantum3/stiffnessplot.xml tutorial/quantum3/stiffnessplot2.xml \
tutorial/qwl/parm1 tutorial/qwl/parm2 tutorial/qwl/parm3 tutorial/qwl/parm4 \
tutorial/dmrg/particle_in_box/parm1 tutorial/dmrg/particle_in_box/parm2 \
tutorial/dmrg/particle_in_box/parm3 \
tutorial/dmrg/particle_in_box/lattices_dmrg.xml \
tutorial/fulldiag/parm1 tutorial/fulldiag/parm2 tutorial/fulldiag/parm3 \
tutorial/fulldiag/parmV15 tutorial/fulldiag/v15-graph.xml tutorial/fulldiag/dd-graph.xml tutorial/fulldiag/model-dspin.xml \
tutorial/sparsediag/parm1 tutorial/sparsediag/parm2 tutorial/sparsediag/parm3 \
tutorial/dmrg/dmrg/parm1 tutorial/dmrg/dmrg/parm2 tutorial/dmrg/dmrg/parm3 \
tutorial/dmrg/dmrg/parm4 tutorial/dmrg/dmrg/parm5

.PHONY : diag dmrg mc qmc
diag dmrg mc qmc :
	$(MAKE) -C $@ all

.PHONY : update-preamble
update-preamble:
	if test -x "$(top_srcdir)/config/update_preamble"; then \
	  find $(top_srcdir) -type f -name '*\.h' -or -name '*\.C' -or -name '*\.hpp' -or -name '*\.h\.in' | xargs $(top_srcdir)/config/update_preamble; \
	fi

.PHONY : test tests example examples
test tests example examples : diag qmc
	$(MAKE) -C test check

.PHONY : release snapshot
release snapshot :
	rm -rf $(distdir)
	$(MAKE) dist
	tar zxf $(distdir).tar.gz
	rm -f $(distdir).tar.gz
	cd $(distdir); rm -rf $(REMOVE_FOR_RELEASE)
	if test "$@" = "snapshot"; then \
	  date=`date +%Y%m%d`; \
	  rm -rf alps-applications-$$date; \
	  mv -f $(distdir) alps-applications-$$date; \
	  tar zcf alps-applications-$$date.tar.gz alps-applications-$$date; \
	  rm -rf alps-applications-$$date; \
	else \
	  tar zcf $(distdir).tar.gz $(distdir); \
	  rm -rf $(distdir); \
	fi
