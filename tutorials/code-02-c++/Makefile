## Uncomment the following line and set the correct path if you haven't set ALPS_ROOT globally
# ALPS_ROOT = /path/to/your/alps/installation

include $(ALPS_ROOT)/share/alps/include.mk

all:	ising

ising:	ising.o
	$(CXX) -o ising ising.o $(LDFLAGS) $(LIBS)

clean:
	$(RM) $(RMFLAGS) ising.o ising

ising.o:	ising.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@
