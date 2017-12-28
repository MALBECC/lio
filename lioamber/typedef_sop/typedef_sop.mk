######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += Sets_smat.f90
INCLUDES += Gets_smat.f90
INCLUDES += Sets_orthog.f90
INCLUDES += Gets_orthog_2m.f90
INCLUDES += Gets_orthog_4m.f90
INCLUDES += Gets_eigens_m.f90
INCLUDES += Gets_eigens_v.f90
INCLUDES += Rebase_bothxx.f90
INCLUDES += Drop_ldvals.f90
INCLUDES += Calc_smat.f90
INCLUDES += Calc_fulldcmp_4m.f90
INCLUDES += Calc_fulldcmp_7m.f90

$(OBJPATH)/typedef_sop.o : $(INCLUDES) typedef_sop.mk
######################################################################
