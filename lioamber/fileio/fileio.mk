######################################################################
included :=
included += restart_commons.f90
included += restart_coef.f90
included += restart_fock.f90
included += restart_rho.f90

$(OBJPATH)/fileio.o : $(included) fileio.mk
######################################################################
