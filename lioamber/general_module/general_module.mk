######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += read_list.f read_list_h.f
INCLUDES += atmorb.f atmorb_h.f
INCLUDES += vector_selection.f
INCLUDES += sdcmp_cholesky.f sdiag_canonical.f

$(obj_path)/general_module.o : $(INCLUDES) general_module.mk
######################################################################
