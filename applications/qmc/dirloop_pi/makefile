include /Users/tamama/Programs/ALPS/share/alps/include.mk

dirloop_pi : main.cpp lowa.cpp lowa.hpp lowa.element.hpp BH_parm.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dirloop_pi lowa.cpp main.cpp $(LDFLAGS) $(LIBS)
dirloop_pi_mpi : main.cpp lowa.cpp lowa.hpp lowa.element.hpp BH_parm.hpp
	$(CXX) $(CPPFLAGS_MPI) $(CXXFLAGS) -o dirloop_pi_mpi lowa.cpp main.cpp $(LDFLAGS_MPI) $(LIBS_MPI)
dirloop_pi_pvm : main.cpp lowa.cpp lowa.hpp lowa.element.hpp BH_parm.hpp
	$(CXX) $(CPPFLAGS_PVM) $(CXXFLAGS) -o dirloop_pi_pvm lowa.cpp main.cpp $(LDFLAGS_PVM) $(LIBS_PVM)
