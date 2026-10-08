CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# The contributor's folder: their .cpp files and manifest.yml. GENERATED gets
# the algorithm.hpp wrapper and selection.hpp that generate.py writes from the
# manifest; it is searched first, so a stray header elsewhere cannot shadow it.
SUBMISSION ?= ../submission
GENERATED  ?= generated
SUB_SRCS    = $(wildcard $(SUBMISSION)/*.cpp)
SUB_HDRS    = $(wildcard $(SUBMISSION)/*.hpp $(SUBMISSION)/*.h)

algorithm_test: main.cpp problems.cpp problem.hpp solver.hpp $(GENERATED)/.stamp $(SUB_SRCS) $(SUB_HDRS)
	$(CXX) $(CXXFLAGS) -I$(GENERATED) -I. -I$(SUBMISSION) main.cpp problems.cpp $(SUB_SRCS) -o $@

$(GENERATED)/.stamp: generate.py $(SUBMISSION)/manifest.yml $(SUB_HDRS)
	python3 generate.py $(SUBMISSION) $(GENERATED)
	@touch $@

.PHONY: clean
clean:
	rm -rf algorithm_test $(GENERATED)
