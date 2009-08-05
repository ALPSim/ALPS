#General ALPS include directory
include ./Makefile.inc
CXXFLAGS += -g -DWERNER -DRUBTSOV -I..

SOURCES=externalsolver.C hirschfyesim.C fouriertransform.C hilberttransformer.C alps_solver.C clustertransformer.C clusterhelper.C 2dsimpson.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

OBJECTS=externalsolver.o hirschfyesim.o fouriertransform.o hilberttransformer.o auxiliaryfunctions.o alps_solver.o clustertransformer.o clusterhelper.o 2dsimpson.o
OBJECTS_DMFT=main.o selfconsistency.o $(OBJECTS)
OBJECTS_HIRSCHFYE=solver_main.o $(OBJECTS)
OBJECTS_INTERPOLATE=interpolate.o $(OBJECTS)
RUBTSOV_DIR=../DMFT-solvers/Rubtsov
WERNER_DIR=../DMFT-solvers/Werner
IGOR_DIR=../../Igor/k_int_mats
OBJECTS_WERNER=$(WERNER_DIR)/impurity.o  $(WERNER_DIR)/fourpoint.o
#OBJECTS_WERNERB1=$(WERNER_DIR)/impurity_beta1.o  $(WERNER_DIR)/fourpoint.o

OBJECTS_RUBTSOV=$(RUBTSOV_DIR)/rubtsov.o $(RUBTSOV_DIR)/rubtsov_io.o $(RUBTSOV_DIR)/rubtsov_solver.o $(RUBTSOV_DIR)/rubtsov_auxiliary.o \
$(RUBTSOV_DIR)/rubtsov_splines.o $(RUBTSOV_DIR)/rubtsov_measurements.o $(RUBTSOV_DIR)/rubtsov_wang_landau.o $(RUBTSOV_DIR)/rubtsov_fastupdate.o \
$(RUBTSOV_DIR)/rubtsov_model.o $(RUBTSOV_DIR)/rubtsov_worm.o $(RUBTSOV_DIR)/rubtsov_fourpoint.o $(RUBTSOV_DIR)/rubtsov_observables.o\
$(RUBTSOV_DIR)/rubtsov_selfenergy.o
OBJECTS_IGOR=$(IGOR_DIR)/a1.o $(IGOR_DIR)/cdln2.o $(IGOR_DIR)/ek22idn.o $(IGOR_DIR)/fexit.o $(IGOR_DIR)/hamilt.o $(IGOR_DIR)/iprint.o $(IGOR_DIR)/applysym.o $(IGOR_DIR)/cekord.o $(IGOR_DIR)/ek2idn.o $(IGOR_DIR)/fopna.o $(IGOR_DIR)/i1mach.o $(IGOR_DIR)/kintegration.o $(IGOR_DIR)/arbek.o $(IGOR_DIR)/cg.o $(IGOR_DIR)/ek3idn.o $(IGOR_DIR)/gfprodctn.o $(IGOR_DIR)/ietol.o $(IGOR_DIR)/cdln1.o $(IGOR_DIR)/creorder.o $(IGOR_DIR)/errmsg.o $(IGOR_DIR)/grnatm.o $(IGOR_DIR)/inverc.o $(IGOR_DIR)/wrrdmtrx.o
OBJECTS_SOLVERS=$(OBJECTS_RUBTSOV) $(OBJECTS_WERNER)

all: dmft dmft_mpi hirschfye interpolate

all-igor: dmft_igor dmft_igor_mpi 

dmft : $(OBJECTS_DMFT) werner rubtsov
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS)
energy : energy.o auxiliaryfunctions.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o energy energy.o auxiliaryfunctions.o $(LDFLAGS) $(LIBS)
dmft_mpi : $(OBJECTS_DMFT) werner rubtsov
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft_mpi $(OBJECTS_DMFT) $(OBJECTS_SOLVERS) $(LDFLAGS_MPI) $(LIBS_MPI)
dmft_igor : $(OBJECTS_DMFT) werner rubtsov igor
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft_igor $(OBJECTS_DMFT) $(OBJECTS_IGOR) $(OBJECTS_SOLVERS) $(LDFLAGS) $(LIBS) -L//sw/lib/gcc4.3/lib/ -lgfortran 
dmft_igor_mpi : $(OBJECTS_DMFT) werner rubtsov igor
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CT_ALGOS) -o dmft_igor_mpi $(OBJECTS_DMFT) $(OBJECTS_IGOR) $(OBJECTS_SOLVERS) $(LDFLAGS_MPI) $(LIBS_MPI) -L//sw/lib/gcc4.3/lib/ -lgfortran 
hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)
interpolate : $(OBJECTS_INTERPOLATE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o interpolate $(OBJECTS_INTERPOLATE) $(LDFLAGS) $(LIBS)
compute_sigma: compute_sigma.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o compute_sigma compute_sigma.cpp $(LDFLAGS) $(LIBS)
hybridization : hybridization.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hybridization hybridization.cpp fouriertransform.o clustertransformer.o clusterhelper.o auxiliaryfunctions.o $(LDFLAGS) $(LIBS)
rubtsov:
	cd $(RUBTSOV_DIR); $(MAKE)
werner:
	cd $(WERNER_DIR); $(MAKE)
igor:
	cd $(IGOR_DIR); $(MAKE)
clean:
	$(RM) $(RMFLAGS) *.o dmft dmft_mpi hirschfye $(RUBTSOV_DIR)/*.o $(WERNER_DIR)/*.o
fileclean:
	rm G_omega_[1-9]* G0_omega_[1-9]* G_tau_[1-9]* G0_tau_[1-9]*
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $<

%.o: %.C
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<
depend:
				makedepend $(SOURCES) $(SOURCES_DMFT)
hirschfyesim.o: hirschfyeaux.h hirschfyesim.h
selfconsistency.o: hilberttransformer.h fouriertransform.h 
fouriertransform.o: fouriertransform.h
main.o: hilberttransformer.h fouriertransform.h
2dsimpson.o: 2dsimpson.C 2dsimpson.h
