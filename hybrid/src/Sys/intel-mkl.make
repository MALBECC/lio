SIESTA_ARCH=intel-mkl
#
# Intel fortran compiler for linux with mkl optimized blas and lapack
#
# Be sure to experiment with different optimization options.
# You have quite a number of combinations to try...
#
FC=ifc 
#
FFLAGS= -w -mp -tpp5 -O3
FFLAGS_DEBUG= -g 
LDFLAGS=-Vaxlib 
COMP_LIBS=
RANLIB=echo
#
NETCDF_LIBS=
NETCDF_INTERFACE=
DEFS_CDF=
#
MPI_INTERFACE=
MPI_INCLUDE=
DEFS_MPI=
#
GUIDE=/opt/intel/mkl/lib/32/libguide.a
LAPACK=/opt/intel/mkl/lib/32/libmkl_lapack.a
BLAS=/opt/intel/mkl/lib/32/libmkl_p3.a
#G2C=/usr/lib/gcc-lib/i386-redhat-linux/2.96/libg2c.a
LIBS=$(LAPACK) $(BLAS) $(G2C) $(GUIDE)  -lpthread 
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








