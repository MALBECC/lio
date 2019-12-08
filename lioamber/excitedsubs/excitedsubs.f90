module excitedsubs
implicit none

contains
#include "ExcProp.f90"
#include "fcaApp.f90"
#include "basis_init.f90"
#include "tda.f90"
#include "initvec.f90"
#include "vecMOtomatMO.f90"
#include "matMOtomatAO.f90"
#include "solve_focks.f90"
#include "MtoIANV.f90"
#include "addInt.f90"
#include "diagonH.f90"
#include "residual.f90"
#include "new_vectors.f90"
#include "norma.f90"
#include "QRfactorization.f90"
#include "OscStr.f90"
#include "TransDipole.f90"
#include "ObtainOsc.f90"
#include "PrintResults.f90"




end module excitedsubs
