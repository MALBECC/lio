################################################################################
INCLUDES :=
INCLUDES += subm_int1.o subm_int1G.o
INCLUDES += subm_int2.o subm_int3lu.o subm_int3mem.o subm_intfld.o subm_intsol.o

OBJECTS += $(INCLUDES)
$(OBJPATH)/faint_cpu.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu.mk
################################################################################
