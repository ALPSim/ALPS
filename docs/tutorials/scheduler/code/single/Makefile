## Copyright (C) 20012 TODO: Licence needet!

## Change ALPS_ROOT if you have installed ALPS in another location
ALPS_ROOT = /opt/alps

include $(ALPS_ROOT)/share/alps/include.mk

all:	ising

ising: ising.o
    $(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ ising.o $(LDFLAGS) $(LIBS)
clean:
    $(RM) $(RMFLAGS) *.o ising
