include $(HOME)/ALPS_INTEL/share/alps/include.mk
#BLAS classes directory
BLAS_CLASSES_DIR=$(HOME)/alps_svn/
#CXX optimization flags
#CXXFLAGS +=  -g -O3 -funroll-loops -ftree-vectorize -march=prescott -W -Wno-sign-compare
#CXX warning flags
#CXXFLAGS += -Wno-deprecated 

CXXFLAGS += -I$(BLAS_CLASSES_DIR)
CXXFLAGS += -DWERNER -DRUBTSOV

