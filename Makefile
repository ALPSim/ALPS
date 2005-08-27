include $(HOME)/ALPS/share/alps/include.mk

SOURCES=hirschfye.C externalsolver.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

OBJECTS=hirschfye.o externalsolver.o
OBJECTS_DMFT=main.o selfconsistency.o $(OBJECTS)
OBJECTS_HIRSCHFYE=solver_main.o $(OBJECTS)

all: dmft hirschfye

dmft : $(OBJECTS_DMFT)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft $(OBJECTS_DMFT) $(LDFLAGS) $(LIBS)

hirschfye : $(OBJECTS_HIRSCHFYE)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye $(OBJECTS_HIRSCHFYE) $(LDFLAGS) $(LIBS)
# DO NOT DELETE

externalsolver.o: externalsolver.h solver.h
hirschfye.o: hirschfye.h solver.h hirschfyesim.h
main.o: selfconsistency.h solver.h hilberttransformer.h hirschfye.h
main.o: externalsolver.h
selfconsistency.o: selfconsistency.h solver.h hilberttransformer.h
solver_main.o: hirschfye.h solver.h
