CHECK =
if BUILD_DIAG_FULLDIAG
CHECK += check
endif

.PHONY : test tests example examples check
test tests example examples : $(CHECK)

OUTXML = parm.task1.out.xml parm.task2.out.xml parm.task3.out.xml parm.task4.out.xml parm.task5.out.xml parm.task6.out.xml parm.task7.out.xml parm.task8.out.xml parm.task9.out.xml parm.task10.out.xml parm.task11.out.xml parm.task12.out.xml parm.task13.out.xml parm.task14.out.xml parm.task15.out.xml parm.task16.out.xml parm.task17.out.xml parm.task18.out.xml parm.task19.out.xml parm.task20.out.xml parm.task21.out.xml parm.task22.out.xml parm.task23.out.xml parm.task24.out.xml

OUTXML_BOSON= boson.task1.out.xml boson.task2.out.xml boson.task3.out.xml boson.task4.out.xml boson.task5.out.xml boson.task6.out.xml boson.task7.out.xml boson.task8.out.xml boson.task9.out.xml boson.task10.out.xml boson.task11.out.xml boson.task12.out.xml boson.task13.out.xml boson.task14.out.xml boson.task15.out.xml boson.task16.out.xml boson.task17.out.xml boson.task18.out.xml boson.task19.out.xml boson.task20.out.xml

OUTXML_FERMION = fermion.task1.out.xml fermion.task2.out.xml fermion.task3.out.xml fermion.task4.out.xml fermion.task5.out.xml fermion.task6.out.xml fermion.task7.out.xml fermion.task8.out.xml fermion.task9.out.xml fermion.task10.out.xml fermion.task11.out.xml fermion.task12.out.xml fermion.task13.out.xml fermion.task14.out.xml fermion.task15.out.xml fermion.task16.out.xml fermion.task17.out.xml fermion.task18.out.xml fermion.task19.out.xml fermion.task20.out.xml fermion.task21.out.xml

OUTXML_FERMION2 = fermion2.task1.out.xml fermion2.task2.out.xml fermion2.task3.out.xml fermion2.task4.out.xml fermion2.task5.out.xml fermion2.task6.out.xml fermion2.task7.out.xml fermion2.task8.out.xml fermion2.task9.out.xml

parm:
	$(LN_S) ../parm .

fermion:
	$(LN_S) ../parm_fermion fermion

boson:
	$(LN_S) ../parm_boson boson

check: check_fermion check_spin check_boson

check_fermion: fermion
	$(DESTDIR)$(bindir)/convert2xml fermion
	../../diag/fulldiag/fulldiag fermion.in.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML_FERMION) > archive.xml
	$(DESTDIR)$(bindir)/xslttransform archive.xml ../energy.xsl > result_fermion.txt

check_boson: boson
	$(DESTDIR)$(bindir)/convert2xml boson
	../../diag/fulldiag/fulldiag boson.in.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML_BOSON) > archive.xml
	$(DESTDIR)$(bindir)/xslttransform archive.xml ../energy.xsl > result_boson.txt

check_spin: parm
	$(DESTDIR)$(bindir)/convert2xml parm
	../../diag/fulldiag/fulldiag parm.in.xml
	$(DESTDIR)$(bindir)/archivecat $(OUTXML) > archive.xml
	$(DESTDIR)$(bindir)/xslttransform archive.xml ../energy.xsl > result.txt
