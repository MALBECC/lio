######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += cumat_init_fin.f90
INCLUDES += cumat_multiply.f90
INCLUDES += cumat_bchange.f90
INCLUDES += cumat_set.f90
INCLUDES += cumat_get.f90
INCLUDES += cumat_misc.f90

$(OBJPATH)/typedef_cumat.o : $(INCLUDES) typedef_cumat.mk
######################################################################
