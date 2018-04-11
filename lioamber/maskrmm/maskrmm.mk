######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += rmmput_dens.f90
INCLUDES += rmmget_dens.f90
INCLUDES += rmmput_fock.f90
INCLUDES += rmmget_fock.f90

INCLUDES += rmmcalc0_init.f90
INCLUDES += rmmcalc2_focknuc.f90
INCLUDES += rmmcalc3_fockele.f90

$(OBJPATH)/maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
