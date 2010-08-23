#General ALPS include directory
include ./Makefile.inc
CXXFLAGS += -g -DWERNER -DRUBTSOV -I..

SOURCES=2dsimpson.C alps_solver.C auxiliaryfunctions.C clusterhelper.C clustertransformer.C externalsolver.C fouriertransform.C hilberttransformer.C\
interpolator.C  main.C selfconsistency.C solver_main.C

OBJECTS_DMFT=2dsimpson.o clusterhelper.o fouriertransform.o selfconsistency.o\
alps_solver.o clustertransformer.o hilberttransformer.o interpolator.o\
auxiliaryfunctions.o externalsolver.o hirschfyesim.o main.o

RUBTSOV_DIR=../DMFT-solvers/Rubtsov
WERNER_DIR=../DMFT-solvers/Werner

IGOR_DIR=../../Igor/k_int_mats
OBJECTS_WERNER=$(WERNER_DIR)/impurity.o  $(WERNER_DIR)/fourpoint.o

OBJECTS_RUBTSOV=$(RUBTSOV_DIR)/rubtsov.o $(RUBTSOV_DIR)/rubtsov_io.o $(RUBTSOV_DIR)/rubtsov_solver.o $(RUBTSOV_DIR)/rubtsov_auxiliary.o \
$(RUBTSOV_DIR)/rubtsov_splines.o $(RUBTSOV_DIR)/rubtsov_measurements.o $(RUBTSOV_DIR)/rubtsov_wang_landau.o $(RUBTSOV_DIR)/rubtsov_fastupdate.o \
$(RUBTSOV_DIR)/rubtsov_model.o $(RUBTSOV_DIR)/rubtsov_worm.o $(RUBTSOV_DIR)/rubtsov_fourpoint.o $(RUBTSOV_DIR)/rubtsov_observables.o\
$(RUBTSOV_DIR)/rubtsov_selfenergy.o
OBJECTS_IGOR=$(IGOR_DIR)/a1.o $(IGOR_DIR)/cdln2.o $(IGOR_DIR)/ek22idn.o $(IGOR_DIR)/fexit.o $(IGOR_DIR)/hamilt.o $(IGOR_DIR)/iprint.o $(IGOR_DIR)/applysym.o $(IGOR_DIR)/cekord.o $(IGOR_DIR)/ek2idn.o $(IGOR_DIR)/fopna.o $(IGOR_DIR)/i1mach.o $(IGOR_DIR)/kintegration.o $(IGOR_DIR)/arbek.o $(IGOR_DIR)/cg.o $(IGOR_DIR)/ek3idn.o $(IGOR_DIR)/gfprodctn.o $(IGOR_DIR)/ietol.o $(IGOR_DIR)/cdln1.o $(IGOR_DIR)/creorder.o $(IGOR_DIR)/errmsg.o $(IGOR_DIR)/grnatm.o $(IGOR_DIR)/inverc.o $(IGOR_DIR)/wrrdmtrx.o
OBJECTS_SOLVERS=$(OBJECTS_RUBTSOV) $(OBJECTS_WERNER)

all: dmft 

all-igor: dmft_igor dmft_igor_mpi 

dmft : $(OBJECTS_DMFT) werner rubtsov
				$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS)
energy : energy.o auxiliaryfunctions.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o energy energy.o auxiliaryfunctions.o $(LDFLAGS) $(LIBS)
dmft_igor : $(OBJECTS_DMFT) werner rubtsov igor
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft_igor $(OBJECTS_DMFT) $(OBJECTS_IGOR) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS) -L//sw/lib/gcc4.3/lib/ -lgfortran 
dmft_igor_mpi : $(OBJECTS_DMFT) werner rubtsov igor
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft_igor_mpi $(OBJECTS_DMFT) $(OBJECTS_IGOR) $(OBJECTS_SOLVERS) $(LDFLAGS_MPI) $(LIBS_MPI) -L//sw/lib/gcc4.3/lib/ -lgfortran 
hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)
compute_sigma: compute_sigma.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o compute_sigma compute_sigma.cpp $(LDFLAGS) $(LIBS)
rubtsov:
	cd $(RUBTSOV_DIR); $(MAKE)
werner:
	cd $(WERNER_DIR); $(MAKE)
igor:
	cd $(IGOR_DIR); $(MAKE)
clean:
	$(RM) $(RMFLAGS) *.o dmft hirschfye $(RUBTSOV_DIR)/*.o $(WERNER_DIR)/*.o
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $<
%.o: %.C
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<
depend:
				makedepend $(SOURCES) $(SOURCES_DMFT)
hirschfyesim.o: hirschfyeaux.h hirschfyesim.h
selfconsistency.o: hilberttransformer.h fouriertransform.h 
main.o: hilberttransformer.h fouriertransform.h
2dsimpson.o: 2dsimpson.C 2dsimpson.h
# DO NOT DELETE
