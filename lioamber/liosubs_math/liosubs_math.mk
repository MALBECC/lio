######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += matdcmp_cholesky.f90
INCLUDES += matdcmp_svd.f90
INCLUDES += matmul3_body.f90
INCLUDES += matmul3_head.f90
INCLUDES += purge_zeros.f90
INCLUDES += transform_gen.f90

$(OBJPATH)/liosubs_math.o : $(INCLUDES) liosubs_math.mk
######################################################################
