#BLAS classes directory
BLAS_CLASSES_DIR=$(HOME)/alps_svn/
#CXX optimization flags
CXXFLAGS +=  -O3 -funroll-loops -ftree-vectorize
#CXX warning flags
CXXFLAGS += -W -Wno-deprecated

CXXFLAGS += -I$(BLAS_CLASSES_DIR)
