SIESTA_ARCH=osfdxml-mpich
#
# The exact name and location of the libraries might vary.
#
FC=mpif90
FFLAGS= -fast -tune host
FFLAGS_DEBUG= -g3 -fast -tune host
RANLIB=echo
LIBS= /usr/local/BLACS/LIB/blacs_MPI-ALPHA-0.a \
      /usr/local/BLACS/LIB/blacsCinit_MPI-ALPHA-0.a \
      /usr/local/BLACS/LIB/blacsF77init_MPI-ALPHA-0.a \
      /usr/local/SCALAPACK/pblas_ALPHA.a \
      /usr/local/SCALAPACK/scalapack_ALPHA.a \
      /usr/local/BLACS/LIB/blacs_MPI-ALPHA-0.a \
      /usr/local/SCALAPACK/pblas_ALPHA.a \
      /usr/local/SCALAPACK/tools_ALPHA.a \
      -ldxml
SYS=bsd
MPI_INTERFACE=libmpi_f90.a
MPI_INCLUDE=/usr/local/include
DEFS=-DMPI
MPILIB=
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

