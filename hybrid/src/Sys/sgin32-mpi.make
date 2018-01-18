SIESTA_ARCH=sgin32-mpi
#
# This file seems to work for SGI systems using the precompiled -n32
# Scalapack and Blacs libraries from netlib, and *some version* (perhaps
# SGI's own?) of MPI.
#
# Note that the Scalapack and Blacs library files must be linked from
# their standard places to the building directory...
#
FC=f90 -n32
#
FFLAGS=  -O3 
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -O0
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/include
DEFS_MPI=-DMPI
#
LIBS= -lcomplib.sgimath libmpiblacsn32.a libscalapackn32.a \
       libmpiblacsn32.a -lmpi 
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $<
#








