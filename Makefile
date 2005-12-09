include $(HOME)/ALPS/share/alps/include.mk

SOURCES=hirschfye.C externalsolver.C hirschfyesim.C fouriertransform.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

OBJECTS=hirschfye.o externalsolver.o hirschfyesim.o fouriertransform.o
OBJECTS_DMFT=main.o selfconsistency.o $(OBJECTS)
OBJECTS_HIRSCHFYE=solver_main.o $(OBJECTS)
OBJECTS_RUBTSOV=rubtsov.o rubtsov_io.o rubtsov_solver.o rubtsov_auxiliary.o rubtsov_splines.o

all: dmft hirschfye

dmft : $(OBJECTS_DMFT)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft $(OBJECTS_DMFT) $(LDFLAGS) $(LIBS)

hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)

rubtsov : $(OBJECTS_RUBTSOV)
	g++ -o rubtsov $(OBJECTS_RUBTSOV) $(LDFLAGS) $(LIBS)
rubtsov_direct : $(OBJECTS_RUBTSOV)
	g++ -o rubtsov_direct $(OBJECTS_RUBTSOV) $(LDFLAGS) $(LIBS)
docs:
	doxygen doxygen.conf
	cd doc/latex; $(MAKE)
	
clean:
	$(RM) $(RMFLAGS) *.o dmft hirschfye
	
# DO NOT DELETE

externalsolver.o: externalsolver.h solver.h
fouriertransform.o: fouriertransform.h hirschfyesim.h solver.h
hirschfye.o: hirschfye.h solver.h hirschfyesim.h
hirschfyesim.o: hirschfyesim.h solver.h
main.o: selfconsistency.h solver.h hilberttransformer.h fouriertransform.h
main.o: hirschfyesim.h hirschfye.h externalsolver.h
selfconsistency.o: selfconsistency.h solver.h hilberttransformer.h
selfconsistency.o: fouriertransform.h hirschfyesim.h
solver_main.o: hirschfye.h solver.h
