#  Copyright (C) Vito Scarola 2013.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)


## Change ALPS_ROOT if you have installed ALPS in another location
ALPS_ROOT = /opt/alps

include $(ALPS_ROOT)/share/alps/include.mk

all:	fqheed

fqheed: fqheed.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ fqheed.o $(LDFLAGS) $(LIBS)  
clean:
	$(RM) $(RMFLAGS) *.o fqheed
