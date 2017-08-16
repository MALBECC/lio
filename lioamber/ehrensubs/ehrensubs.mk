######################################################################
included :=
included += ehrendyn.f90
included += ehrenstep.f90
included += ehrensetup.f90
included += setim.f90
included += calc_forceDS.f90
included += calc_forceDS_dss.f90
included += calc_forceDS_dds.f90

included += ehrenrst.f90
included += ehren_setfld.f90
included += ehren_magnus.f90
included += ehren_verlet.f90
included += ehren_masses.f90
included += calc_kenergy.f90

included += RMMcalc0_Init.f90
included += RMMcalc1_Overlap.f90
included += RMMcalc2_FockMao.f90
included += RMMcalc3_FockMao.f90
included += RMMcalc4_FockMao.f90

$(OBJPATH)/ehrensubs.o : $(included) ehrensubs.mk
######################################################################
