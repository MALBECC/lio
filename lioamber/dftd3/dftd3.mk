included :=
#included += dftd3_3.f90
#included += dftd3_2.f90
#included += dftd3_read_params.f90
#included += dftd3_main.f90

#dftd3_main.f90: dftd3_read_params.o dftd3_2.o dftd3_3.o
dftd3.mod: dftd3_data.mod

$(OBJPATH)/dftd3.o : $(included) dftd3.mk