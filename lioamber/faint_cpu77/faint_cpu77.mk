################################################################################
INCLUDES :=
INCLUDES += subm_int2G.o
INCLUDES += subm_int3G.o

OBJECTS += $(INCLUDES)
$(OBJPATH)/subm_int3G.o  : $(OBJPATH)/subm_int2G.o
$(OBJPATH)/faint_cpu77.o : $(INCLUDES:%.o=$(OBJPATH)/%.o) faint_cpu77.mk

################################################################################
# PRIVATE FLAGS FOR OPTIMIZATION:
#
# The following section assigns to specific OBJECTS certain flags, such as
# which level of optimization corresponds to each object. Intermediate flags
# are used to store the private information and then these option is added
# to the general compilation flags.

#$(OBJPATH)/int1G.o    : private_flag := -O1
#$(OBJPATH)/int3G.o    : private_flag := -O1
#$(OBJPATH)/intfld.o   : private_flag := -O1
#$(OBJPATH)/intsol.o   : private_flag := -O1
#$(OBJPATH)/intsolG.o  : private_flag := -O1

#$(OBJPATH)/int3lu.o   : private_flag := -O3 $(OPTIMP)

#$(OBJPATH)/int1.o     : private_flag := -O3 $(OPTIMI)
#$(OBJPATH)/int2.o     : private_flag := -O3 $(OPTIMI)
#$(OBJPATH)/int2G.o    : private_flag := -O3 $(OPTIMI)
#$(OBJPATH)/intSG.o    : private_flag := -O3 $(OPTIMI)
#$(OBJPATH)/int3mem.o  : private_flag := -O3 $(OPTIMI)

#FFLAGS += $(private_flag)

################################################################################
