SIESTA_ARCH=t3e
#
# Makefile include file for Cray's T3E
# Alberto Garcia <wdpgaara@lg.ehu.es>, Aug 2, 1999
#
# Compiler invocation. Note '-em' to produce module files and
# '-dp' to disable double precision
#
FC=f90 -em -dp
#
# Whatever needed, except -dp and -em... 
#
FFLAGS=  -O scalar2,pipeline2,aggress
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -Rabc -ei
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
DEFS_MPI=-DMPI
#
LIBS= -lsci -lmpi 
#LIBS= -lsci -lmpi -lapp
#
SYS=t3e
RANLIB=echo
DEFS= $(DEFS_CDF) $(DEFS_MPI) -DCRAY
#
# Actual compilation recipes for siesta code.
# Specify "-p ." to let the compiler know that
# the modules are in the current directory.
# It could work without it...
#
.F.o:
	$(FC) -c $(FFLAGS) -p . $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) -p .  $<
.F90.o:
	$(FC) -c $(FFLAGS) -p .  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS) -p .  $<
#

