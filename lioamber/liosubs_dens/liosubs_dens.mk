######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += simple_guess.f90
INCLUDES += builds_densmat.f90
INCLUDES += messup_densmat.f90

$(OBJPATH)/liosubs_dens : $(INCLUDES) liosubs_dens.mk
######################################################################
