included :=

included += TSHmain.f90
included += obtain_wavefunction.f90
included += obtain_sigma.f90
included += fstsh_init.f90
included += fstsh_forces.f90
included += print_tsh.f90
included += do_electronic_interpolation.f90
included += match_CIS.f90

$(OBJPATH)/fstshsubs.o: $(included) fstshsubs.mk

