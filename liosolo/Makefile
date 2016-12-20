LIBS=
FLAGS = -fopenmp 
ifeq ($(intel),1)
FLAGS += -fpp
else
FLAGS += -cpp
endif

ifeq ($(xeon_phi),1)
  FLAGS += -mmic
  LIBS= -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/mic -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/mic/
endif
######################################################################
ifeq ($(intel),1)
  FC   = ifort
  LIBS+= -module ../lioamber/obj -I../lioamber/obj 
else
  FC   = gfortran
  LIBS+= -I../lioamber/obj -J../lioamber/obj
endif

LIBS += -L../g2g -lg2g
LIBS += -L../lioamber -llio-g2g

ifeq ($(magma),1)
        DEFINE   += -Dmagma
        LIBS     += -L$(MAGMAROOT)/lib -lmagma
endif

#FORTRAN_WRAPER=
ifeq ($(cublas),1)
#LIBS += -lstdc++ -L$(CUDA_HOME)/lib64/ -lcudart -lcublas
#FORTRAN_WRAPER=../lioamber/obj/fortran.o
endif
#FLAGS +=$(FORTRAN_WRAPER)


######################################################################
all: liosolo liomd.x

liosolo: liosolo.f90 ../lioamber/liblio-g2g.so
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

liomd.x: liomd.f90 ../lioamber/liblio-g2g.so
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf liosolo liomd.x
######################################################################
