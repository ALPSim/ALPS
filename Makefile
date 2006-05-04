include $(HOME)/ALPS/share/alps/include.mk
CXXFLAGS +=-W -funroll-loops  -g
SOURCES=externalsolver.C hirschfyesim.C fouriertransform.C alps_solver.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

OBJECTS=externalsolver.o hirschfyesim.o fouriertransform.o auxiliaryfunctions.o alps_solver.o
OBJECTS_DMFT=main.o selfconsistency.o $(OBJECTS)
OBJECTS_HIRSCHFYE=solver_main.o $(OBJECTS)
RUBTSOV_DIR=../DMFT-solvers/Rubtsov
WERNER_DIR=../DMFT-solvers/Werner
OBJECTS_WERNER=$(WERNER_DIR)/impurity.o 

OBJECTS_RUBTSOV=$(RUBTSOV_DIR)/rubtsov.o $(RUBTSOV_DIR)/rubtsov_io.o $(RUBTSOV_DIR)/rubtsov_solver.o $(RUBTSOV_DIR)/rubtsov_auxiliary.o $(RUBTSOV_DIR)/rubtsov_splines.o $(RUBTSOV_DIR)/rubtsov_measurements.o $(RUBTSOV_DIR)/rubtsov_wang_landau.o $(RUBTSOV_DIR)/rubtsov_fastupdate.o $(RUBTSOV_DIR)/rubtsov_model.o
OBJECTS_SOLVERS=$(OBJECTS_RUBTSOV) $(OBJECTS_WERNER)

all: dmft dmft_mpi hirschfye

dmft : $(OBJECTS_DMFT) werner rubtsov
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS)

dmft_mpi : $(OBJECTS_DMFT) rubtsov
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft_mpi $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS_MPI) $(LIBS_MPI)

hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)
rubtsov:
	cd ../DMFT-solvers/Rubtsov; $(MAKE)
werner:
	cd ../DMFT-solvers/Werner; $(MAKE)
four:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o four four.C $(LDFLAGS) fouriertransform.o auxiliaryfunctions.o $(LIBS)
docs:
	doxygen doxygen.conf
	cd doc/latex; $(MAKE)
test: dmft
	cd test && make all	
clean:
	$(RM) $(RMFLAGS) *.o dmft hirschfye $(RUBTSOV_DIR)/*.o $(WERNER_DIR)/*.o
fileclean:
	rm G_omega_[1-9]* G0_omega_[1-9]* G_tau_[1-9]* G0_tau_[1-9]*
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

