######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += cumfx.f90
INCLUDES += cumpxt.f90 cumpxt_h.f90
INCLUDES += cumpxt_r.f90
INCLUDES += cumxtf.f90
INCLUDES += cumpx.f90  cumpx_h.f90
INCLUDES += cumpx_r.f90
INCLUDES += cumxp.f90 cumxp_h.f90
INCLUDES += cumxp_r.f90
INCLUDES += cumxtp.f90 cumxtp_h.f90
INCLUDES += cu_fock_commuts.f90

$(OBJPATH)/cublasmath.o : $(INCLUDES) cublasmath.mk

######################################################################
