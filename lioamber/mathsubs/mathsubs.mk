######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += commutator.f90 commutator_h.f90
INCLUDES += basechange_gemm.f90 basechange_gemm_h.f90

$(OBJPATH)/mathsubs.o : $(INCLUDES) mathsubs.mk
######################################################################
