include $(HOME)/ALPS/share/alps/include.mk
#BLAS classes directory
BLAS_CLASSES_DIR=$(HOME)/alps_svn/applications/dmft
#CXX optimization flags

#CXX warning flags
CXXFLAGS += -Wall -Wno-deprecated 

CXXFLAGS += -I$(BLAS_CLASSES_DIR)
CXXFLAGS += -DWERNER 

