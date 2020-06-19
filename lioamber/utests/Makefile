######################################################################
# MAKEFILE FOR TESTS
######################################################################

FC = gfortran
TESTS = test-properties.x

all: $(TESTS)

ljs:
	$(FC) test-ljswitch.f90 -o ljs.x ../lennard_jones/ljs_mm_interface.f90 -cpp -Wall -g
 
%.x: %.f90
	$(FC) $< -o $@

.PHONY: clean
clean:
	rm -rf *.x *.out *.mod

######################################################################
