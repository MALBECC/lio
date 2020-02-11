######################################################################
# MAKEFILE FOR TESTS
######################################################################
FC := gfortran

tests :=
tests += matdiag_ut

all : $(tests:%=%.x)

and_run : $(tests:%=%.out)
	@echo ""
	@echo "------RESULTS------"
	@more all.out
	@rm all.out
	@echo ""


%.x : %.f90 Makefile
	$(FC) $< -o $@

%.out : %.x Makefile
	@./$< > $@
	@cat $@ >> all.out
	@echo >> all.out
	@rm $@

.PHONY: clean
clean:
	rm -rf $(tests:%=%.x)

######################################################################
