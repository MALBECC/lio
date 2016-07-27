######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += rmmput_dens_all.f90  rmmput_dens_h.f90
INCLUDES += rmmget_dens_all.f90  rmmget_dens_h.f90
INCLUDES += rmmput_fock_all.f90  rmmput_fock_h.f90
INCLUDES += rmmget_fock_all.f90  rmmget_fock_h.f90
INCLUDES += rmmcalc_fockmao.f90  rmmcalc_overlap.f90


$(obj_path)/maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
