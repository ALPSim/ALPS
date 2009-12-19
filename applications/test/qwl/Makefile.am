CHECK =
if BUILD_QMC_QWL
CHECK += check
endif

.PHONY : test tests example examples check
test tests example examples : $(CHECK)

dist_noinst_DATA = parm

check: parm
	$(DESTDIR)$(bindir)/convert2xml parm
	../../qmc/qwl/qwl --Tmin 2 parm.in.xml
	../../qmc/qwl/qwl_evaluate parm.task1.out.xml
	../../qmc/qwl/qwl_evaluate parm.task2.out.xml
	$(DESTDIR)$(bindir)/plot2text parm.task1.plot.energy.xml > result.txt
	$(DESTDIR)$(bindir)/plot2html parm.task2.plot.entropy.xml > result.html

