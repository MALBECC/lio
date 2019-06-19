######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += cumat_bchange.f90
INCLUDES += cumat_blas_subs.f90
INCLUDES += cumat_exchange.f90
INCLUDES += cumat_get.f90
INCLUDES += cumat_init_fin.f90
INCLUDES += cumat_misc.f90
INCLUDES += cumat_multiply.f90
INCLUDES += cumat_set.f90
INCLUDES += cumat_update.f90

$(OBJPATH)/typedef_cumat.o : $(INCLUDES) typedef_cumat.mk
######################################################################
