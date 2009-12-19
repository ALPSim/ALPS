CHECK =
if BUILD_QMC_LOOPER
CHECK += check
endif

.PHONY : test tests example examples check
test tests example examples : $(CHECK)

OUTXML = parm.task1.out.xml parm.task2.out.xml parm.task3.out.xml parm.task4.out.xml parm.task5.out.xml parm.task6.out.xml parm.task7.out.xml parm.task8.out.xml parm.task9.out.xml parm.task10.out.xml parm.task11.out.xml parm.task12.out.xml parm.task13.out.xml parm.task14.out.xml parm.task15.out.xml parm.task16.out.xml parm.task17.out.xml parm.task18.out.xml parm.task19.out.xml parm.task20.out.xml parm.task21.out.xml parm.task22.out.xml parm.task23.out.xml parm.task24.out.xml

parm:
	$(LN_S) ../parm .

check: parm
	$(DESTDIR)$(bindir)/parameter2xml parm
	../../qmc/looper/loop parm.in.xml
	../../qmc/looper/loop --evaluate parm.out.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML) > archive.xml
	$(DESTDIR)$(bindir)/xslttransform archive.xml ../energy.xsl > result.txt
