################################################################################
# -lname tells the linker to look for external procedures in libname library.
#
# -Lpathname tells the linker to look for libraries in pathname.
#
# -Wl,-rpath='pathname' saves pathname inside the executable so it will look
# for libraries there at runtime
#
LIBS  =
LIBS += -Wl,-rpath='$$ORIGIN/../g2g' -L../g2g -lg2g
LIBS += -Wl,-rpath='$$ORIGIN/../lioamber' -L../lioamber -llio-g2g
ifneq ($(cuda),0)
LIBS += -L$(CUDA_HOME)/lib64 -lcudart -lcuda
endif


FLAGS = -fopenmp


################################################################################
ifeq ($(intel),1)
  FC    = ifort
  FLAGS+= -fpp
  LIBS += -module ../lioamber/obj -I../lioamber/obj
else ifeq ($(intel),2)
  FC    = ifort
  FLAGS+= -fpp
  LIBS += -module ../lioamber/obj -I../lioamber/obj
else
  FC    = gfortran
  FLAGS+= -cpp
  LIBS += -I../lioamber/obj -J../lioamber/obj
endif

ifeq ($(xeon_phi),1)
  FLAGS+= -mmic
  LIBS += -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/mic
  LIBS += -Qoption,link,-rpath=/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/mic/
endif


ifeq ($(magma),1)
  DEFINE += -Dmagma
  LIBS   += -L$(MAGMAROOT)/lib -lmagma
endif


################################################################################
all: liosolo liomd.x

liosolo: liosolo.f90 ../lioamber/liblio-g2g.so
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

liomd.x: liomd.f90 ../lioamber/liblio-g2g.so
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf liosolo liomd.x

################################################################################
