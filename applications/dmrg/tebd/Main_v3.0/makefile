#Makefile for the OSTEBD-ALPS interface
#----------------------------------------------
HDF5DIR=/usr/local/hdf5
OFLAGS = -O3 -ffree-line-length-none -fopenmp#Set the optimization flags
LIBFLAGS = -llapack -lblas #Link the LAPACK and BLAS libraries
FC=gfortran#/opt/local/bin/gfortran-mp-4.4
HDF5Flags = -L$(HDF5DIR)/lib -I$(HDF5DIR)/include -lhdf5 -lhdf5_fortran -Wl,-rpath,$(HDF5DIR)/lib -lz
MODDIR=../Core_v3.0/


#-----------------------------------
#List objects
#-----------------------------------
OBJS = $(MODDIR)GlobalData.o $(MODDIR)LinearOps.o $(MODDIR)HamiOps.o $(MODDIR)StateOps.o  $(MODDIR)LocalOps.o $(MODDIR)ObOps.o $(MODDIR)PropOps.o  $(MODDIR)Hdf5Interface.o

#-----------------------------------
#Object compilation rules
#-----------------------------------

moddir_obj :
	cd $(MODDIR);  make all

#Main Programs	
tebd.o : tebd.f90
		$(FC) $(OFLAGS) -c tebd.f90 -I$(MODDIR) $(HDF5Flags)


#-----------------------------------
#Primary rules
#-----------------------------------

tebd :  moddir_obj tebd.o
	$(FC) $(OFLAGS)  -o tebd tebd.o $(OBJS) $(LIBFLAGS) -I$(MODDIR) $(HDF5Flags)


clean: #removes all old *.o and *.mod files
	rm *.o; rm tebd;
	cd $(MODDIR) ;	rm *.o ;	rm *.mod 


	
