# arch.make for ifort compiler

SIESTA_ARCH=ifc-libs
FC=ifort
FFLAGS= -mp1 -ip -O3 
LDFLAGS= -static-libcxa
RANLIB=ranlib
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FREE_F90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $(FREE_F90_CPP) $<
