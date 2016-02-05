######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += set_masses.f90 writegeom.f90 nuclear_verlet.f90

$(obj_path)/mathsubs.o : $(INCLUDES) liosubs.mk
######################################################################
