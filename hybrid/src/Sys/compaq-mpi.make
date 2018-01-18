# Makefile include file for parallel Compaq machine at Grenoble

SIESTA_ARCH=compaq-mpi

FC=f90 

FFLAGS=  -O2
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -Rabc -ei

NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=

MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
DEFS_MPI=-DMPI

LIBS = -L/usr/local/lib/scalapack -lscalapack -lpblas -ltools -lblacsF77 \
        -lblacs -lblacsF77 \
        -lblacs -ldxml -lfmpi -lmpi -lelan

SYS=bsd
RANLIB=echo
DEFS= $(DEFS_CDF) $(DEFS_MPI) 

# Actual compilation recipes for siesta code.

.F.o:
        $(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
        $(FC) -c $(FFLAGS)   $<
.F90.o:
        $(FC) -c $(FFLAGS)   $(DEFS) $<
.f90.o:
        $(FC) -c $(FFLAGS)   $<

