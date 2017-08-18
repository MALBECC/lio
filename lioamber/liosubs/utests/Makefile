######################################################################
# MAKEFILE FOR TESTS
######################################################################
FC     := gfortran
FFLAGS := -cpp

tests :=
tests += nuclear_verlet_t01.x
tests += find_free_unit_ut.x

all : $(tests)

%.x : %.f90 Makefile
	$(FC) $(FFLAGS) $< -o $@
	./$@
	rm -rf $@ *.mod

.PHONY: clean
clean:
	rm -rf $(tests) *.mod

######################################################################
