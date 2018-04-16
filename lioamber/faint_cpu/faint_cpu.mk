################################################################################
INCLUDES :=
INCLUDES += subm_int1.o

OBJECTS += $(INCLUDES)
#$(OBJPATH)/faint_cpu.o : $(INCLUDES) faint_cpu.mk
$(OBJPATH)/faint_cpu.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu.mk
################################################################################
