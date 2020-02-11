######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += gaussbell.f90 gaussbell_h.f90
INCLUDES += commutator.f90 commutator_h.f90
INCLUDES += basechange.f90 basechange_h.f90
INCLUDES += basechange_gemm.f90 basechange_gemm_h.f90

$(OBJPATH)/mathsubs.o : $(INCLUDES) mathsubs.mk
######################################################################
