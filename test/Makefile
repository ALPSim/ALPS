all: difftest test
	echo "tests done. exiting."
difftest: difftest.cpp
	g++ -o difftest difftest.cpp

test: difftest
	cd rubtsov_pm && make
	cd hilbert && make 
	cd no_half_filling && make
	cd hirschfye && make
#	cd rubtsov_wl && make test
