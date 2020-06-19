INCLUDES := lj_switch_data.o
OBJECTS += $(INCLUDES)

$(OBJPATH)/lj_switch.o: lj_switch_data.o ljs_dEdQ.f90 ljs_init_end.f90 ljs_mm_interface.f90 ljs_fock_terms.f90 ljs_input.f90
$(OBJPATH)/lj_switch.o: $(INCLUDES) ljswitch.mk