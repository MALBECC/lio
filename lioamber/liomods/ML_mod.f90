#include "datatypes/datatypes.fh"
module ML_mod 
    implicit none 

    LIODBLE, allocatable    :: Tmat_vec(:)              ! Vectorized form of the kinetic energy matrix
    LIODBLE                 :: kinE = 0.D0              ! Kinetic energy

   ! Variable to control output level when generating density fitting data for ML models.
   integer :: df_verbosity   ! 0 = none, 1 = energies+coeffs, 2 = +derivatives

   ! Variable to control whether the kinetic energy is computed or nto
   logical :: compute_kinE

end module ML_mod
