######################################################################
# INTERNAL DEPENDENCIES
######################################################################
internal_files := liosubs.mk
internal_files += liosubs.f90
internal_files += find_free_unit.f90
internal_files += catch_iostat.f90
internal_files += set_masses.f90
internal_files += nuclear_verlet.f90
internal_files += write_energy.f90
internal_files += write_geom.f90
$(obj_path)/liosubs.o : $(internal_files)

######################################################################
# EXTERNAL DEPENDENCIES
######################################################################
external_users := SCF.o
$(external_users:%.o=$(obj_path)/%.o) : $(obj_path)/liosubs.mod

######################################################################
