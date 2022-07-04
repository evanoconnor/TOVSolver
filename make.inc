F90=gfortran
F90FLAGS= -O3 -g -ffixed-line-length-none
MODINC="-I ./"
LDFLAGS= -O3 -g -ffixed-line-length-none

HDF5DIR=/Users/eoco/opt/hdf5-1.13.0-gfortran
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5_fortran -lhdf5 -lz
HDF5INCS=-I$(HDF5DIR)/include

