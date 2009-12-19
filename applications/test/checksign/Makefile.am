CHECK =
if BUILD_QMC_CHECKSIGN
CHECK += check
endif

.PHONY : test tests example examples check
test tests example examples : $(CHECK)

check:
	../../qmc/checksign/checksign ../parm
