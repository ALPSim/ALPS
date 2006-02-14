

all:
	difftest
	test
difftest: difftest.cpp
	g++ -o difftest difftest.cpp

test: difftest
	cd rubtsov_pm && make test
#	cd rubtsov_wl && make test
