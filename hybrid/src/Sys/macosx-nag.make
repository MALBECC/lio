SIESTA_ARCH=macosx-nag-cdf
#
#  NAG Fortran compiler on MacOS X
#  
# Issues addressed: 
#
#      - Case-insensitive filesystem
#      - .F and .F90 suffixes not recognized by f95 on Mac (possibly related
#        to the first point. Explicit pre-processing.
#      - CPUtime now called in the f95 fashion (file nag.f)
#      - -dcfuns option to allow non-standard intrinsics such as dconjg
#      - -dusty to allow for some deprecated features in blas and lapack.
# 
# Other issues:
#
#      ONE HAS TO BE EXTRA CAREFUL WHEN DELETING FILES, DUE TO THE
#      CASE-INSENSITIVENESS OF THE FILE SYSTEM ON THE MAC.
#
FC=f95
FC_ASIS=$(FC)
RANLIB=ranlib
#
FFLAGS= -g -dcfuns -dusty
FFLAGS_DEBUG= 
LDFLAGS=
COMP_LIBS=
#
NETCDF_LIBS=            # -L/opt/lib -lnetcdf
NETCDF_INTERFACE=       # libnetcdf_f90.a
DEFS_CDF= -D__NAG__     # -DCDF
FFLAGS_NETCDF=-mismatch_all -w=unused          # Relax module checks...
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
LIBS=  $(MPI_LIBS)  $(NETCDF_LIBS) #-framework veclib  (maybe for G5's)
COMP_LIBS=linalg.a
SYS=nag
DEFS= $(DEFS_CDF) $(DEFS_MPI)
#
CPP=/usr/local/lib/NAGWare/fpp -P
#
%.o : %.mod
#
.F.o:
	$(CPP) -fixed $(DEFS) $*.F > aux_$*.f
	$(FC) -c $(FFLAGS)  $(DEFS) aux_$*.f 
	@rm -f aux_$*.f
	@mv aux_$*.o $*.o
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(CPP) -free $(DEFS) $*.F90 > aux_$*.f90
	$(FC) -c $(FFLAGS) aux_$*.f90
	@rm -f aux_$*.f90
	@mv aux_$*.o $*.o
.f90.o:
	$(FC) -c $(FFLAGS)   $<
#









