SIESTA_ARCH=macosx-xlf
#
# IBM XL Fortran Compiler on MacOS X  (beta version)
#
# Issues addressed: 
#
#      - Case-insensitive filesystem
# Other issues:
#
#      ONE HAS TO BE EXTRA CAREFUL WHEN DELETING FILES, DUE TO THE
#      CASE-INSENSITIVENESS OF THE FILE SYSTEM ON THE MAC.
#
FC=xlf
FC_ASIS=$(FC)
FREE_F90=-qsuffix=f=f90 -qfree=f90
FREE_F90_CPP=-qsuffix=cpp=F90 -qfree=f90
#
RANLIB=ranlib
#
FFLAGS=-O
FFLAGS_DEBUG= -g -O0
LDFLAGS=
COMP_LIBS=
#
NETCDF_LIBS=          #-L/opt/lib -lnetcdf
NETCDF_INTERFACE=     # libnetcdf_f90.a
DEFS_CDF=             # -WF,-DCDF
FFLAGS_NETCDF= -qsuffix=f=f90:cpp=F90 
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS=  $(MPI_LIBS)  $(NETCDF_LIBS) 
COMP_LIBS=linalg.a
SYS=xlf
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
%.o : %.mod
#
# Beta compiler crashes on these...
#
atom.o:
	$(FC) -c $(FFLAGS_DEBUG)  atom.f
siesta.o:
	$(FC) -c  $(FFLAGS_DEBUG)  $(DEFS)  siesta.F
#
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(FC) -c $(FREE_F90_CPP) $(FFLAGS)  $(DEFS) $<
.f90.o:
	$(FC) -c $(FREE_F90) $(FFLAGS)   $<
#









