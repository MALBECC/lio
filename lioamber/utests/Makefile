######################################################################
# MAKEFILE FOR TESTS
######################################################################

LIBS=

FLAGS = -fopenmp
ifeq ($(intel),1)
    FLAGS += -fpp
else
    FLAGS += -cpp
endif

# Libraries and compiler.
ifeq ($(intel),1)
  FC   = ifort
  LIBS+= -module ../obj -I../obj
else
  FC   = gfortran
  LIBS+= -I../obj -J../obj
endif

LIBS += -L../../g2g -lg2g
LIBS += -L../ -llio-g2g

tests :=
tests += test-properties.x

all : $(tests)

%.x : %.f Makefile
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

%.x : %.F Makefile
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

%.x : %.f90 Makefile
	$(FC) $< -o $@ $(FLAGS) $(LIBS)

%.x : %.F90 Makefile
	$(FC) $< -o $@ $(FLAGS) $(LIBS)


.PHONY: clean
clean:
	rm -rf $(tests)

######################################################################
