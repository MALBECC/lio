SIESTA_ARCH=sun-nolibs
#
# Serial compilation without the need of any installed libraries.
# SUN version
#
FC=f90
FC_ASIS=$(FC)
RANLIB=echo
#
FFLAGS= -O
FFLAGS_DEBUG= -g -O0
LIBS=  
SYS=bsd
DEFS=
MPILIB=
COMP_LIBS=linalg.a
MPI_INCLUDE=/usr/local/include
CPP=/usr/local/bin/cpp -P 
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(FFLAGS) $*.f90
	@rm -f $*.f90
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
# %.o : %.mod
