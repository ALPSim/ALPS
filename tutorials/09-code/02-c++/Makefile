#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

## Change ALPS_ROOT if you have installed ALPS in another location
ALPS_ROOT?= /opt/alps

include $(ALPS_ROOT)/share/alps/include.mk

all:	ising

ising: ising-skeleton.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ ising-skeleton.o $(LDFLAGS) $(LIBS)
clean:
	$(RM) $(RMFLAGS) *.o ising
