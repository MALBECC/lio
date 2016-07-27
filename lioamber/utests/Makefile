######################################################################
# MAKEFILE FOR TESTS
######################################################################
FC := gfortran

tests :=
tests += vvfterm_test01.x
tests += lowdinpop_test01.x

all : $(tests)

%.x : %.f Makefile
	$(FC) $< -o $@

%.x : %.F Makefile
	$(FC) $< -o $@

%.x : %.f90 Makefile
	$(FC) $< -o $@

%.x : %.F90 Makefile
	$(FC) $< -o $@

.PHONY: clean
clean:
	rm -rf $(tests)

######################################################################
