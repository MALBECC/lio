######################################################################
SRCDIRS += ehrenfest

# INTERNAL DEPENDENCIES
included :=
included += ehrendyn.f90
included += ehrensetup.f90
included += setim.f90
included += calc_forceDS.f90
included += calc_forceDS_dss.f90
included += calc_forceDS_dds.f90

included += ehrenrst.f90
included += ehren_magnus.f90
included += ehren_verlet_e.f90
included += ehren_masses.f90
included += calc_kenergy.f90

included += RMMcalc0_Init.f90
included += RMMcalc1_Overlap.f90
included += RMMcalc2_FockMao.f90

$(OBJPATH)/ehrenfest.o : $(included) ehrenfest.mk
######################################################################
