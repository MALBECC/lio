######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += builds_densmat.f90
INCLUDES += messup_densmat.f90
INCLUDES += starting_guess.f90
INCLUDES += standard_coefs.f90

$(OBJPATH)/liosubs_dens.o : $(INCLUDES) liosubs_dens.mk
######################################################################
