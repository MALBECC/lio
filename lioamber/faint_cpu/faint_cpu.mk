################################################################################
INCLUDES :=
INCLUDES += subm_int1.o subm_int2.o subm_int3lu.o subm_int3mem.o subm_intsol.o
INCLUDES += subm_int1G.o subm_intSG.o subm_intsolG.o
INCLUDES += subm_intfld.o

OBJECTS += $(INCLUDES)
$(OBJPATH)/faint_cpu.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu.mk
################################################################################
