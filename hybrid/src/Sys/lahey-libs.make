SIESTA_ARCH=lahey-libs
#
FC=lf95
FC_ASIS=$(FC)
#
FFLAGS= -O3 --warn --quiet --tpp7 --ntrace
FFLAGS_DEBUG= -g -O0  --chk --trace
LDFLAGS= --staticlink
RANLIB=echo
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
LIBS=
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








