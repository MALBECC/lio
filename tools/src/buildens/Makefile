.DEFAULT_GOAL := buildens

ifeq ($(intel),1)
   FC = ifort
else ifeq ($(intel),2)
   FC = ifort
else
   FC = gfortran
endif


buildens:
	$(FC) -o buildens.x buildens.f90
	cp buildens.x ../../bdens

clean:
	rm -f buildens.x ../../bdens
