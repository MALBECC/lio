######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += cuconmut_r.f   
INCLUDES += cumatmulnanoc.f cumatmulnanoc_h.f 
INCLUDES += cumfx.f  
INCLUDES += cumpxt.f cumpxt_h.f
INCLUDES += cumpxt_r.f
INCLUDES += cumxtf.f 
INCLUDES += cuconmut.f cuconmut_h.f     
INCLUDES += cumagnusfac.f cumagnusfac_h.f 
INCLUDES += cumatmul_r.f     
INCLUDES += cumpx.f  cumpx_h.f 
INCLUDES += cumpx_r.f   
INCLUDES += cumxp.f cumxp_h.f
INCLUDES += cumxp_r.f  
INCLUDES += cumxtp.f cumxtp_h.f 
INCLUDES += cupredictor.f cupredictor_h.f


$(obj_path)/cublasmath.o : $(INCLUDES) cublasmath.mk
######################################################################
