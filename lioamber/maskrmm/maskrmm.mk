######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += spunpack.f sprepack.f
INCLUDES += rhomess.f  rhomess_h.f
INCLUDES += rhofix.f   rhofix_h.f

$(obj_path)/maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
