######################################################################
included :=
included += restart_commons.f90
included += restart_coef.f90
included += restart_fock.f90
included += restart_rho.f90
included += restart_td.f90
included += output_init.f90
included += output_scf.f90

restart_coef.o restart_rho.o restart_fock.o: restart_commons.o
restart_td.o: restart_rho.o restart_commons.o

$(OBJPATH)/fileio.o : $(included) fileio.mk
######################################################################
