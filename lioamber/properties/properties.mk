INCLUDES := properties_data.o
OBJECTS += $(INCLUDES)

$(OBJPATH)/properties.o: properties_data.o becke.f90 fukui.f90 lowdin.f90 misc.f90 \
                         mulliken.f90 write_population.f90 dipole.f90 init_fin.f90
$(OBJPATH)/properties.o: $(INCLUDES:%.o=$(OBJPATH)/%.o) properties.mk
