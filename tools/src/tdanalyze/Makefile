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
	cp tdanalyze.x ../../tdanalyze
	cp average_spectra.x ../../average_spectra

clean:
	rm -f average_spectra.x ../../average_spectra
	rm -f tdanalyze.x ../../tdanalyze
