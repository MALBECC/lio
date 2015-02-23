######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += gaussbell.f gaussbell_h.f
INCLUDES += commutator.f commutator_h.f
INCLUDES += basechange.f basechange_h.f

$(obj_path)/mathsubs.o : $(INCLUDES) mathsubs.mk
######################################################################
