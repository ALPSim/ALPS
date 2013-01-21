#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

## Change ALPS_ROOT if you have installed ALPS in another location
ALPS_ROOT = /opt/alps

include $(ALPS_ROOT)/share/alps/include.mk

all: ising_c

ising_c: ising.o export2py.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -bundle -o $@.so ising.o export2py.o $(LDFLAGS) $(LIBS)
run:
	PYTHONPATH=$PYTHONPATH:. $(ALPS_ROOT)/bin/alpspython main.py -T 10 result.h5
clean:
	$(RM) $(RMFLAGS) *.o ising_c.so result.h5
