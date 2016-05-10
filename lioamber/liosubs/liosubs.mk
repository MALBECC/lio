######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += find_free_unit.f90
INCLUDES += catch_iostat.f90
INCLUDES += set_masses.f90
INCLUDES += nuclear_verlet.f90

$(obj_path)/liosubs.o : $(INCLUDES) liosubs.mk
######################################################################
