################################################################################
INCLUDES :=
INCLUDES += subm_int1.o subm_int1G.o

OBJECTS += $(INCLUDES)
$(OBJPATH)/faint_cpu.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu.mk
################################################################################
