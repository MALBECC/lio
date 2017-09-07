######################################################################
included :=
included += ehrendyn_main.f90
included += ehrendyn_init.f90
included += ehrendyn_step.f90
included += ehrendyn_prep.f90

included += calc_gintmat.f90
included += calc_forceDS.f90
included += calc_forceDS_dss.f90
included += calc_forceDS_dds.f90
included += calc_Dmat_cholesky.f90
included += calc_Dmat_h.f90

included += ehrenaux_rstx.f90
included += ehrenaux_setfld.f90
included += ehrenaux_magnus.f90
included += ehrenaux_verlet.f90
included += ehrenaux_masses.f90
included += ehrenaux_calckyn.f90
included += ehrenaux_writedip.f90
included += ehrenaux_updatevel.f90

included += RMMcalc0_Init.f90
included += RMMcalc1_Overlap.f90
included += RMMcalc2_FockMao.f90
included += RMMcalc3_FockMao.f90
included += RMMcalc4_FockMao.f90

$(OBJPATH)/ehrensubs.o : $(included) ehrensubs.mk
######################################################################
