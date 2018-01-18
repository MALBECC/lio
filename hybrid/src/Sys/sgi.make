SIESTA_ARCH=sgi
#
FC=f90 -n32
RANLIB=echo
#
FFLAGS=  -O3 
#FFLAGS= -O scalar2,pipeline2,aggress -eA
FFLAGS_DEBUG= -g -O0
#
LIBS= -lcomplib.sgimath 
SYS=bsd
#
# Location of mpif.h include file
#
MPILIB=
MPI_INCLUDE=/usr/local/include
#
# Definition to trigger MPI conditional compilation
#
DEFS=
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
