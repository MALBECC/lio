######################################################################
# MAKEFILE FOR TESTS
######################################################################

FC = gfortran
TESTS = test-properties.x test-ftermbias.x

all: $(TESTS)

%.x: %.f90
	$(FC) $< -o $@

.PHONY: clean
clean:
	rm -rf *.x *.out

######################################################################
