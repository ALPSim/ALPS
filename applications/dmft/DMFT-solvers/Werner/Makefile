DMFT_DIR= ../../DMFT/
include $(DMFT_DIR)/Makefile.inc
CXXFLAGS+= -I../..

all: impurity.o fourpoint.o

impurity.o: moves.h update.h impurity.h impurity.C 
fourpoint.o: moves.h update.h impurity.h fourpoint.C 
clean:
	$(RM) ti_files *.o *~ 
depend:
				makedepend -I../../ impurity.C fourpoint.C
