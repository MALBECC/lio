included := dftd3.o

dftd3.mod: dftd3_data.mod

$(OBJPATH)/dftd3.o : $(included) dftd3.mk