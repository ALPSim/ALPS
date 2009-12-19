CHECK =
if BUILD_QMC_SSE
CHECK += check
endif

.PHONY : test tests example examples check
test tests example examples : $(CHECK)

OUTXML = parm.task1.out.xml parm.task2.out.xml parm.task3.out.xml parm.task4.out.xml parm.task5.out.xml parm.task6.out.xml parm.task7.out.xml parm.task8.out.xml parm.task9.out.xml parm.task10.out.xml parm.task11.out.xml parm.task12.out.xml parm.task13.out.xml parm.task14.out.xml parm.task15.out.xml parm.task16.out.xml parm.task17.out.xml parm.task18.out.xml parm.task19.out.xml parm.task20.out.xml parm.task21.out.xml parm.task22.out.xml parm.task23.out.xml parm.task24.out.xml

OUTXML_BOSON = boson.task1.out.xml boson.task2.out.xml boson.task3.out.xml boson.task4.out.xml boson.task5.out.xml boson.task6.out.xml boson.task7.out.xml boson.task8.out.xml boson.task9.out.xml boson.task10.out.xml boson.task11.out.xml boson.task12.out.xml boson.task13.out.xml boson.task14.out.xml boson.task15.out.xml boson.task16.out.xml boson.task17.out.xml boson.task18.out.xml boson.task19.out.xml boson.task20.out.xml

parm:
	$(LN_S) ../parm .

boson:
	$(LN_S) ../parm_boson boson

check: check_spin check_boson

check_spin: parm
	$(DESTDIR)$(bindir)/convert2xml parm
	../../qmc/sse/dirloop_sse --Tmin 2 parm.in.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML) > archive.xml
	$(DESTDIR)$(bindir)/xslttransform archive.xml ../energy.xsl > result.txt

check_boson: boson
	$(DESTDIR)$(bindir)/convert2xml boson
	../../qmc/sse/dirloop_sse --Tmin 2 boson.in.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML_BOSON) > archive_boson.xml
	$(DESTDIR)$(bindir)/xslttransform archive_boson.xml ../energy.xsl > result_boson.txt
