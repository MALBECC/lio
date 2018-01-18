SIESTA_ARCH=sgi64-mpi_fermat
#
# This file seems to work for SGI systems using 64-bit compiled
# Scalapack and Blacs libraries from netlib, and *some version* (perhaps
# SGI's own?) of MPI.
#
# Note that the Scalapack and Blacs library files must be linked from
# their standard places to the building directory... 
#
# Note that the locations of the libraries are configured for the
# Fermat system of the CSAR service at Manchester. Details must
# changed for other machines accordingly.
#
FC=f90 -64
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
LIBS= -lscalapack -lblacs -lpblas -ltools \
      /usr/local/unsupported/numerical/lib64/blacsCinit_MPI-O2K_64-0.a \
      /usr/local/unsupported/numerical/lib64/blacs_MPI-O2K_64-0.a \
      -lscs -lcomplib.sgimath -lmpi 
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








