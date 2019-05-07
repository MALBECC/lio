######################################################################
included :=
included += converger_subs.f90
included += converger_commons.f90
included += converger_ls.f90

converger_commons.o: converger_ls.o converger_diis.o
converger_subs.o   : converger_commons.o

$(OBJPATH)/converger_subs.o : $(included) conver.mk
######################################################################
