SIESTA_ARCH=pgf90
#
FC=pgf90
FC_ASIS=$(FC)
#
FFLAGS=  -fast
FFLAGS_DEBUG= -g -O0
RANLIB=echo
LDFLAGS=
COMP_LIBS=dc_lapack.a
#
NETCDF_LIBS=            # /usr/local/netcdf-3.5/lib/pgi/libnetcdf.a
NETCDF_INTERFACE=       # libnetcdf_f90.a
DEFS_CDF=               # -DCDF
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS= -L/usr/local/lib -llapack -lblas -lg2c $(NETCDF_LIBS)
SYS=bsd
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
#
# Important (at least for V5.0-1 of the pgf90 compiler...)
# Compile atom.f without optimization.
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG) atom.f
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
