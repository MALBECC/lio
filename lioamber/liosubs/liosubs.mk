###############################################################################
# INTERNAL DEPENDENCIES
###############################################################################
internal_files := liosubs.mk
internal_files += liosubs.f90

internal_files += catch_error.f90
internal_files += catch_iostat.f90
internal_files += find_free_unit.f90
internal_files += safeio_open.f90
internal_files += safeio_rewind.f90
internal_files += write_energy.f90
internal_files += write_geom.f90

internal_files += set_masses.f90
internal_files += nuclear_verlet.f90

$(OBJPATH)/liosubs.o : $(internal_files)

###############################################################################
