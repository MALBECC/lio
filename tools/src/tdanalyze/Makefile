.DEFAULT_GOAL := tdanalyze

ifeq ($(intel),1)
   FC = ifort
else ifeq ($(intel),2)
   FC = ifort
else
   FC = gfortran
endif


tdanalyze:
	$(FC) -o tdanalyze.x tdanalyze.f90
	$(FC) -o average_spectra.x average_spectra.f90

clean:
	rm -rf average_spectra
	rm -rf tdanalyze
