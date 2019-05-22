######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += cumat_init_fin.f90
INCLUDES += cumat_multiply.f90


$(OBJPATH)/typedef_cumat.o : $(INCLUDES) typedef_cumat.mk
######################################################################
