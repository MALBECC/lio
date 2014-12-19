######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += spunpack.f sprepack.f
INCLUDES += rhomess.f  rhomess_h.f
INCLUDES += rhofix.f   rhofix_h.f

maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
