INCLUDES :=
OBJECTS += $(INCLUDES)

$(OBJPATH)/dftd3.o: dftd3_2.f90  dftd3_3.f90  dftd3_c6c8.f90  dftd3.f90  dftd3_main.f90  dftd3_read_params.f90  param_c6.f90  param_r0.f90
$(OBJPATH)/dftd3.o: $(INCLUDES:%.o=$(OBJPATH)/%.o) dftd3.mk