

## Change ALPS_ROOT if you have installed ALPS in another location
ALPS_ROOT = /opt/alps

include $(ALPS_ROOT)/share/alps/include.mk

all:	fqheed

fqheed: fqheed.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ fqheed.o $(LDFLAGS) $(LIBS)  
clean:
	$(RM) $(RMFLAGS) *.o fqheed
