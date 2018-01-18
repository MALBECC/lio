SIESTA_ARCH=osfdxml
#
FC=f90
FFLAGS= -fast -tune host
FFLAGS_DEBUG= -g
RANLIB=echo
COMP_LIBS=
#
NETCDF_LIBS=-L/usr/local/netcdf-3.5/lib -lnetcdf
NETCDF_INTERFACE=libnetcdf_f90.a
DEFS_CDF=-DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= -ldxml $(NETCDF_LIBS)
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
CPP=/bin/cpp -P
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








