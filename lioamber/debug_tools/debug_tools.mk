######################################################################
included :=
included += debug_tools.f90
included += Check_posvel.f90
included += Print_matrix.f90

$(OBJPATH)/debug_tools.o : $(included) debug_tools.mk
######################################################################
