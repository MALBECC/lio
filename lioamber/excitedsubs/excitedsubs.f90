#include "../datatypes/datatypes.fh"
module excitedsubs
implicit none

contains
#include "ExcProp.f90"
#include "fcaApp.f90"
#include "basis_exc.f90"
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
#include "OscStr.f90"
#include "TransDipole.f90"
#include "ObtainOsc.f90"
#include "PrintResults.f90"
#include "calc2eFITT.f90"
#include "RelaxedDensity.f90"
#include "GenerateDensities.f90"
#include "Zvector.f90"
#include "ChangeBasisF.f90"
#include "Rcalculate.f90"
#include "PCG_solve.f90"
#include "forcesexc.f90"
#include "WSgradcalc.f90"
#include "HVgradcalc.f90"
#include "COUgradcalc.f90"
#include "tsh_probabilities.f90"
#include "intSG_Exc.f90"
#include "reduced_space.f90"
#include "truncated_MOs.f90"
#include "energy_specific.f90"
#include "getcubegen_excited.f90"




end module excitedsubs
