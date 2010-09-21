#makefile for Core_v3.0 subdirectory
HDF5DIR=/usr/local/hdf5
OFLAGS = -O3 -fopenmp#Set the optimization flags
LIBFLAGS = -llapack -lblas #Link the LAPACK and BLAS libraries
HDF5Flags = -L$(HDF5DIR)/lib -I$(HDF5DIR)/include -lhdf5 -lhdf5_fortran -Wl,-rpath,$(HDF5DIR)/lib -lz
FC=gfortran#/opt/local/bin/gfortran-mp-4.4

#-----------------------------------
#List objects
#-----------------------------------
OBJS = GlobalData.o LinearOps.o HamiOps.o StateOps.o  LocalOps.o ObOps.o PropOps.o Hdf5Interface.o

#-----------------------------------
#Object compilation rules
#-----------------------------------

%.o : %.f90
	$(FC) $(OFLAGS) -c $< $(HDF5Flags)

all : $(OBJS)
