################################################################################
# -lname tells the linker to look for external procedures in libname library.
#
# -Lpathname tells the linker to look for libraries in pathname.
#
# -Wl,-rpath='pathname' saves pathname inside the executable so it will look
# for libraries there at runtime
#

LIBS_NATIVE =
LIBS_NATIVE += -Wl,-rpath='$$ORIGIN/../g2g' -L../g2g -lg2g
LIBS_NATIVE += -Wl,-rpath='$$ORIGIN/../lioamber' -L../lioamber -llio-g2g

LIBS_LIBXC = 
ifneq ($(libxc),0)
# Add liblio-g2g-libxc.so library for libxc functionality
LIBS_LIBXC += -Wl,-rpath='$$ORIGIN/../g2g' -L../g2g -lg2g_libxc
LIBS_LIBXC += -Wl,-rpath='$$ORIGIN/../lioamber' -L../lioamber -llio-g2g-libxc
endif


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
#all: liosolo liomd.x
all: liosolo 

#liosolo: liosolo.f90 ../lioamber/liblio-g2g.so
#	$(FC) $< -o $@ $(FLAGS) $(LIBS)

# Link together the main executable with both implementations
liosolo: liosolo.o 
	$(FC) $^ -o $@ $(FLAGS) -L. -lliosolo_libxc -lliosolo_native $(LIBS) $(LIBS_NATIVE) $(LIBS_LIBXC) 

# Compile individual object files
liosolo.o: liosolo.f90 libliosolo_native.so libliosolo_libxc.so
	$(FC) $(FLAGS) -c $< -o $@ $(LIBS)

libliosolo_native.so: liosolo_native.f90 ../lioamber/liblio-g2g.so
	$(FC) $(FLAGS) -shared -fPIC $< -o $@ $(LIBS_NATIVE) $(LIBS)

libliosolo_libxc.so: liosolo_libxc.f90 ../lioamber/liblio-g2g.so
	$(FC) $(FLAGS) -shared -fPIC -DUSE_LIBXC=1 $< -o $@ $(LIBS_LIBXC) $(LIBS) 

liomd.x: liomd.f90 ../lioamber/liblio-g2g.so
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

clean:
	rm -rf liosolo liomd.x libliosolo_libxc.so libliosolo_native.so *.o *.mod

################################################################################
