######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += gaussbell.f gaussbell_h.f
INCLUDES += commutate.f commutate_h.f
INCLUDES += transbase.f transbase_h.f

mathsubs.o : $(INCLUDES) mathsubs.mk
######################################################################
