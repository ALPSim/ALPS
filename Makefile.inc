include /opt/alps/share/alps/include.mk
#BLAS classes directory
BLAS_CLASSES_DIR=$(HOME)/alps_svn
#CXX optimization flags
#CXXFLAGS +=  -g -O3 -funroll-loops -ftree-vectorize -march=native -W -Wno-sign-compare
#CXXFLAGS +=  -g -O3 -funroll-loops -ftree-vectorize -march=native -W -Wno-sign-compare
#CXXFLAGS+=-march=native
#CXX warning flags
#CXXFLAGS += -Wno-deprecated 

CXXFLAGS += -I$(BLAS_CLASSES_DIR) -I$(DMFT_ROOT_DIR)
CXXFLAGS += -DWERNER -DRUBTSOV -Wno-sign-compare -Wno-deprecated

