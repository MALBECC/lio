################################################################################
INCLUDES :=
INCLUDES += int1.f int1G.f
INCLUDES += int2.f int2G.f
INCLUDES += int3lu.f int3mem.f int3G.f
INCLUDES += intSG.f intfld.f
INCLUDES += intsol.f intsolG.f intsolGs.f

$(OBJPATH)/faint_cpu77.o : $(INCLUDES) faint_cpu77.mk
################################################################################
