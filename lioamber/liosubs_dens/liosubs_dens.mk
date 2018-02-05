######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += starting_guess.f90
INCLUDES += builds_densmat.f90
INCLUDES += messup_densmat.f90

$(OBJPATH)/liosubs_dens.o : $(INCLUDES) liosubs_dens.mk
######################################################################
