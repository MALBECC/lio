included :=

included += ExcProp.f90
included += fcaApp.f90
included += basis_exc.f90
included += tda.f90
included += initvec.f90
included += vecMOtomatMO.f90
included += matMOtomatAO.f90
included += solve_focks.f90
included += MtoIANV.f90
included += addInt.f90
included += diagonH.f90
included += residual.f90
included += new_vectors.f90
included += norma.f90
included += QRfactorization.f90
included += OscStr.f90
included += TransDipole.f90
included += ObtainOsc.f90
included += PrintResults.f90
included += calc2eFITT.f90



$(OBJPATH)/excitedsubs.o: $(included) excitedsubs.mk
