include $(HOME)/ALPS/share/alps/include.mk

SOURCES=hirschfye.C externalsolver.C
SOURCES_DMFT= main.C selfconsistency.C $(SOURCES) 

all: dmft hirschfye

dmft : $(SOURCES_DMFT)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o dmft $(SOURCES_DMFT) $(LDFLAGS) $(LIBS)

hirschfye : solver_main.C $(SOURCES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o hirschfye solver_main.C $(SOURCES) $(LDFLAGS) $(LIBS)
