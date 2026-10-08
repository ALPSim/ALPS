CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

algorithm_test: main.cpp problems.cpp problem.hpp solver.hpp algorithm.hpp
	$(CXX) $(CXXFLAGS) main.cpp problems.cpp -o algorithm_test

.PHONY: clean
clean:
	rm -f algorithm_test
