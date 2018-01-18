SIESTA_ARCH=darwin
#
FC=f90
FC_ASIS=$(FC)
#
FFLAGS=  -O
FFLAGS_DEBUG= -g
RANLIB=echo "... do we need to ranlib this? :"
LDFLAGS=
COMP_LIBS=linalg.a
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= -lblas -lU77 $(NETCDF_LIBS)
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
.F.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#
