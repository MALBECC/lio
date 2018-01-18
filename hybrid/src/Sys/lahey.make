SIESTA_ARCH=lahey
#
FC=lf95
FC_ASIS=$(FC)
RANLIB=echo
#
FFLAGS= -O --warn --quiet --tpp --ntrace
FFLAGS_DEBUG= -g -O0  --chk --trace
LDFLAGS=--staticlink
COMP_LIBS=
#
NETCDF_LIBS=-L/usr/local/netcdf-3.5/lib/lahey -lnetcdf
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= -L/usr/local/lib/lahey -llapack -lblas $(NETCDF_LIBS)
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








