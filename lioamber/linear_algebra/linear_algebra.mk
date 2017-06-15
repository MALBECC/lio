######################################################################
# INTERNAL DEPENDENCIES
######################################################################
internal_files := linear_algebra.mk
internal_files += linear_algebra.f90
internal_files += matrix_diagon_d.f90
internal_files += matrix_diagon_dsyevd.f90  matrix_diagon_dsyevr.f90
internal_files += matmult_interface.f90
internal_files += matmult_procedures.f90 matmult_body.f90
$(OBJPATH)/linear_algebra.o : $(internal_files)

######################################################################
# EXTERNAL DEPENDENCIES
######################################################################
object_users := SCF.o
$(object_users:%.o=$(OBJPATH)/%.o) : $(OBJPATH)/linear_algebra.mod

######################################################################
