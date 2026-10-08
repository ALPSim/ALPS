CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# The contributor's folder: algorithm.hpp plus any .cpp files it calls into.
SUBMISSION ?= ../submission
SUB_SRCS    = $(wildcard $(SUBMISSION)/*.cpp)
SUB_HDRS    = $(wildcard $(SUBMISSION)/*.hpp $(SUBMISSION)/*.h)

algorithm_test: main.cpp problems.cpp problem.hpp solver.hpp $(SUB_SRCS) $(SUB_HDRS)
	@test -f $(SUBMISSION)/algorithm.hpp || \
	  { echo "error: $(SUBMISSION)/algorithm.hpp is missing" >&2; exit 1; }
	$(CXX) $(CXXFLAGS) -I. -I$(SUBMISSION) main.cpp problems.cpp $(SUB_SRCS) -o algorithm_test

.PHONY: clean
clean:
	rm -f algorithm_test
