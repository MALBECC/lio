# arch.make for ifort compiler without library linking

GAUSSDFT_ARCH=ifc-nolibs
FC=ifort
FFLAGS= -mp1 -ip -O3
RANLIB=ranlib
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FREE_F90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $(FREE_F90_CPP) $<
