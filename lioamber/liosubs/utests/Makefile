######################################################################
# MAKEFILE FOR TESTS
######################################################################
FC := gfortran

tests :=
tests += nuclear_verlet_t01.x
tests += find_free_unit_t01.x

all : $(tests)

%.x : %.f Makefile
	$(FC) $< -o $@

%.x : %.f90 Makefile
	$(FC) $< -o $@

.PHONY: clean
clean:
	rm -rf $(tests)

######################################################################
