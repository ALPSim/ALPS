#General ALPS include directory
include ./Makefile.inc
CXXFLAGS += -DWERNER 

SOURCES=externalsolver.C hirschfyesim.C fouriertransform.C hilberttransformer.C alps_solver.C 2dsimpson.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

OBJECTS=externalsolver.o hirschfyesim.o fouriertransform.o hilberttransformer.o auxiliaryfunctions.o alps_solver.o 2dsimpson.o
OBJECTS_DMFT=main.o selfconsistency.o $(OBJECTS)
OBJECTS_HIRSCHFYE=solver_main.o $(OBJECTS)
WERNER_DIR=../DMFT-solvers/Werner
OBJECTS_WERNER=$(WERNER_DIR)/impurity.o  $(WERNER_DIR)/fourpoint.o

OBJECTS_SOLVERS=$(OBJECTS_WERNER)

all: dmft dmft_mpi hirschfye

dmft : $(OBJECTS_DMFT) werner 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS)
dmft_mpi : $(OBJECTS_DMFT) werner 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft_mpi $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS_MPI) $(LIBS_MPI)
hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)
werner:
	cd $(WERNER_DIR); $(MAKE)
clean:
	$(RM) $(RMFLAGS) *.o dmft dmft_mpi hirschfye $(RUBTSOV_DIR)/*.o $(WERNER_DIR)/*.o
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $<
%.o: %.C
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<
depend:
				makedepend $(SOURCES_DMFT) $(SOURCES_WERNER)
