###############################################################################
# INTERNAL DEPENDENCIES
###############################################################################
internal_files := liosubs.mk
internal_files += liosubs.f90

internal_files += line_search.f90
internal_files += catch_error.f90
internal_files += check_vecsize.header.f90 check_vecsize.proced.f90
internal_files += check_matsize.header.f90 check_matsize.proced.f90
internal_files += atmvec_to_orbvec.header.f90 atmvec_to_orbvec.proced.f90
internal_files += read_list.header.f90 read_list.proced.f90

internal_files += gaussian_shaper.f90
internal_files += set_masses.f90
internal_files += nuclear_verlet.f90

internal_files += catch_iostat.f90
internal_files += find_free_unit.f90
internal_files += safeio_open.f90
internal_files += safeio_rewind.f90
internal_files += write_energy.f90
internal_files += write_geom.f90

$(OBJPATH)/liosubs.o : $(internal_files)

###############################################################################
